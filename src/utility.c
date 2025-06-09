// Contains implementations of utility.h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <utility.h>

void VectorPushBack(Vector* self, void* data)
{
    // resize if needed
    if (self->count == self->limit)
    {
        self->limit = self->limit * 2;  // Doubles the size of the array when it needs more capacity
        self->data = realloc(self->data, sizeof(void*) * self->limit);  // void* is the variable type
    }

    // add data
    self->data[self->count] = data;     // count is index
    self->count ++;     // increment count after
}

void VectorFastRemove(Vector* self, int index)
{
    if (index > -1 && index < self->count)
    {
        self->data[index] = self->data[self->count - 1]; // Swap index with last item
        self->data[self->count - 1] = NULL; // Last item is made null
        self->count --;
    }
    
}

void VectorFree(Vector* self)
{
    if (self->data) // There is data
    {
        free(self->data);
        self->data = NULL;
    }
}

void VectorInit(Vector* vector)
{
    vector->limit = LIMIT;
    vector->count = 0;
    vector->push_back = VectorPushBack;
    vector->fast_remove = VectorFastRemove;
    vector->free = VectorFree;
    vector->data = malloc(sizeof(void*)*vector->limit);
}

int comp(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

double median(double array[])
{
    int size = sizeof(array) / sizeof(array[0]);
    qsort(array, size, sizeof(double), comp);

    if (size % 2 != 0)
    {
        return array[size/2];
    }

    else
    {
        double ans = (array[(size - 1)/2] + array[size/2]) / 2.0;
        return ans;
    } 
    
}