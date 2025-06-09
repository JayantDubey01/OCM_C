#ifndef UTILITY_H
#define UTILITY_H

#define LIMIT 100

typedef struct Vector 
{
    void **data;    // ** = Pointer to any data type 
    int limit;      // Tracking capacity
    int count;      // Tracking size
    int dim;        // Initialize 1D, 2D, or 3D dynamically allocated memory
    void (*push_back)(struct Vector*, void*);   // Function pointer for adding
    void (*fast_remove)(struct Vector*, int);   // Function pointer for removing
    void (*free)(struct Vector*);               // Function pointer for freeing

} Vector;

void VectorInit(Vector* vector);

double* median(double array[]);

#endif