#include <math.h>

typedef struct {
    float x;
    float y;
    float z;
} Vector;

float dot(Vector *a, Vector *b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

float len(Vector *a) {
    return sqrt(dot(a, a));
}