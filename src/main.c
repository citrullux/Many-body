#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "types.h"

// сила, действующая на частицу p1 со стороны частицы p2
Vector force(Vector *r1, Vector *r2) {
    Vector r;
    r.x = r1->x - r2->x;
    r.y = r1->y - r2->y;
    r.z = r1->z - r2->z;
    float l2 = dot(&r, &r);

    Vector f = {0};
    // оборванный потенциал Леннарда-Джонса (6-12)
    if (l2 < 2.5 * 2.5) {
        float a = 12 * pow(l2, -7) - 6 * pow(l2, -4);
        f.x = r.x * a;
        f.y = r.y * a;
        f.z = r.z * a;
    }
    return f;
}

int main(int argc, char *argv[]) {
    int n = 1000;
    float all_time = 1;
    float dt = 1e-3;
    int rank, size;
    FILE *fd;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* create a type for type Vector */
    const int    nitems=3;
    int          blocklengths[3] = {1,1,1};
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype MPI_VECTOR;
    MPI_Aint     offsets[3];

    offsets[0] = offsetof(Vector, x);
    offsets[1] = offsetof(Vector, y);
    offsets[2] = offsetof(Vector, z);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_VECTOR);
    MPI_Type_commit(&MPI_VECTOR);


    Vector *rs = (Vector *)malloc(n * sizeof(Vector));
    
    int n1 = n / size;
    Vector *vs1 = malloc(n1 * sizeof(Vector));
    Vector *forces1 = malloc(n1 * sizeof(Vector));

    if (rank == 0) {
        Vector *forces = (Vector *)malloc(n * sizeof(Vector));
        Vector *vs = (Vector *)malloc(n * sizeof(Vector));
        srand(time(NULL));
        // начальные условия
        puts("init");
        int l = (int)pow(n, 1 / 3.) + 1;
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < l; ++j) {
                for (int k = 0; k < l; ++k) {
                    int index = i * l * l + j * l + k;
                    if (index < n) {
                        rs[index].x = i - l / 2;
                        rs[index].y = j - l / 2;
                        rs[index].z = k - l / 2;

                        vs[index].x = 2 * (float)rand() / RAND_MAX - 1;
                        vs[index].y = 2 * (float)rand() / RAND_MAX - 1;
                        vs[index].z = 2 * (float)rand() / RAND_MAX - 1;
                    } else { break; }
                }
            }
        }
        // предварительный расчёт скоростей для момента dt/2
        puts("v(dt/2)");
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                Vector f = force(rs + i, rs + j);
                
                forces[i].x += f.x;
                forces[i].y += f.y;
                forces[i].z += f.z;
                
                forces[j].x -= f.x;
                forces[j].y -= f.y;
                forces[j].z -= f.z;
            }
        }
        for (int i = 0; i < n; ++i) {
            vs[i].x += forces[i].x * dt / 2;
            vs[i].y += forces[i].y * dt / 2;
            vs[i].z += forces[i].z * dt / 2;
        }
        fd = fopen("data.out", "wb");
        fwrite(rs, sizeof(Vector), n, fd);
        
        for (int i = 0; i < n1; ++i)
        {
            vs1[i] = vs[i];
        }
        puts("start send");
        for (int i = 1; i < size; ++i) {
            MPI_Send(vs + i * n1, n1, MPI_VECTOR, i, 2, MPI_COMM_WORLD);
        }
        puts("end send");
        free(vs);
        free(forces);
    } else {
        MPI_Status status;
        MPI_Recv(vs1, n1, MPI_VECTOR, 0, 2, MPI_COMM_WORLD, &status);
    }

    // симуляция
    if (rank == 0) {
        puts("simulation");
    }
    for (float t = 0; t < all_time; t += dt) {
        // отправляем положения частиц всем процессам
        if (rank == 0) {
            printf("%f of %f\n", t, all_time);
            for (int i = 1; i < size; ++i) {
                MPI_Send(rs, n, MPI_VECTOR, i, 1, MPI_COMM_WORLD);
            }    
        } else {
            // получаем из главного процесса положения атомов
            MPI_Status status;
            MPI_Recv(rs, n, MPI_VECTOR, 0, 1, MPI_COMM_WORLD, &status);
        }
        
        // перерасчёт положений
        for (int i = 0; i < n1; ++i) {
            rs[rank * n1 + i].x += vs1[i].x * dt;
            rs[rank * n1 + i].y += vs1[i].y * dt;
            rs[rank * n1 + i].z += vs1[i].z * dt;
        }
        // перерасчёт сил для частиц, за которые отвечает этот процесс
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != rank * n1 + i) {
                    Vector f = force(rs + rank * n1 + i, rs + j);
                    forces1[i].x += f.x;
                    forces1[i].y += f.y;
                    forces1[i].z += f.z;
                }
            }
        }
        // перерасчёт скоростей
        for (int i = 0; i < n1; ++i) {
            vs1[i].x += forces1[i].x * dt;
            vs1[i].y += forces1[i].y * dt;
            vs1[i].z += forces1[i].z * dt;
        }
        if (rank != 0) {
            MPI_Send(rs + rank * n1, n1, MPI_VECTOR, 0, 1, MPI_COMM_WORLD);
        }
        if (rank == 0) {
            // получаем данные от процессов
            for (int i = 1; i < size; ++i)
            {
                MPI_Status status;
                MPI_Recv(rs + i * n1, n1, MPI_VECTOR, i, 1, MPI_COMM_WORLD, &status);
            }
            
            // дампим в файл
            fwrite(rs, sizeof(Vector), n, fd);
        }
    }

    if (rank == 0) {
        fclose(fd);
    }
    free(forces1);
    free(rs);
    free(vs1);
    MPI_Finalize();
    return 0;
}