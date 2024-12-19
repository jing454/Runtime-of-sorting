// FILE: select.cxx
// An interactive test program for the selectionsort function

#include <algorithm>  // Provides swap
#include <cstdlib>    // Provides EXIT_SUCCESS, size_t
#include <iostream>   // Provides cout and cin
#include <time.h>     //provides clock_t, clock, CLOCKS_PER_SEC
using namespace std;

void insertionsort(int data[], size_t n);                    //header for insertionsort

void selectionsort(int data[], size_t n);                    //header for selectionsort

void quicksort(int data[], size_t n);                        //header for quicksort 

void partition(int data[], size_t n, size_t& pivot_index);          //header for quicksort part 2

void mergesort(int data[], size_t n);                               //mergesort

void merge(int data[], size_t n1, size_t n2);                       //mergesort part 2

void heapsort(int data[], size_t n);                               //heapsort

size_t parent(size_t k);                                           //heap

size_t left_child(size_t k);                                       //heap
                          
size_t right_child(size_t k);                                      //heap

void make_heap(int data[], size_t n);                              //heap

void reheapify_down(int data[], size_t n);                         //reheapification

int main() {                      //main 

	const size_t n = 500000;                                 //array size
	int ARRAY_SIZE = n;                                       //array size
	size_t number_of_element = n;
	int* data = new int[ARRAY_SIZE];                                    //allocate array on heap
	srand(time(NULL));                               //ensure rand() produce random values
	for (int i = 0; i < number_of_element; i++) {      //inputs random data 1-ARRAY_SIZE into every index of the array
		data[i] = rand() % ARRAY_SIZE + 1;
	}

	clock_t start_insertionsort_unsort = clock();                   //starts clock for unsorted (insertion)

	insertionsort(data, number_of_element);                   //call insortionsort

	clock_t end_insertionsort_unsort = clock();                  //ends clock for unsorted (insertion)

	double clock_insertionsort_unsort = double(end_insertionsort_unsort - start_insertionsort_unsort) / CLOCKS_PER_SEC;     //finds the amount of seconds that it took

	cout << "insertionsort of unsorted took " << clock_insertionsort_unsort << " seconds" << endl;             //prints the time it look to sort an unsorted array (seconds)

	clock_t start_insertionsort_sort = clock();                                          //start clock for sorted

	insertionsort(data, number_of_element);                          //call insortionsort     

	clock_t end_insertionsort_sort = clock();                               //ends clock for sorted                     

	double clock_insertionsort_sort = double(end_insertionsort_sort - start_insertionsort_sort) / CLOCKS_PER_SEC;            //finds the amount of seconds that it took

	cout << "insertionsort of sorted took " << clock_insertionsort_sort << " seconds" << endl;          //prints the time it took to sort an unsorted array (seconds)

	delete[] data;                      //deallocate 

	data = new int[ARRAY_SIZE];              //re-allocate 
	
	for (int i = 0; i < number_of_element; i++) {      //inputs random data 1-ARRAY_SIZE into every index of the array
		data[i] = rand() % ARRAY_SIZE + 1;
	}
	clock_t start_selectionsort_unsort = clock();            //timer for selection sort starts

	selectionsort(data, number_of_element);                   //call selectionsort

	clock_t end_selectionsort_unsort = clock();               //timer for selection sort stops

	double clock_selectionsort_unsort = double(end_selectionsort_unsort - start_selectionsort_unsort) / CLOCKS_PER_SEC;     //finds the amount of seconds that it took

	cout << "selectionsort of unsorted took " << clock_selectionsort_unsort << " seconds" << endl;             //prints the time it look to sort an unsorted array (seconds)

	clock_t start_selectionsort_sort = clock();                    //starts clock for sorted

	selectionsort(data, number_of_element);                     //call selectionsort     

	clock_t end_selectionsort_sort = clock();                 //starts clock for sorted

	double clock_selectionsort_sort = double(end_selectionsort_sort - start_selectionsort_sort) / CLOCKS_PER_SEC;            //finds the amount of seconds that it took

	cout << "selectionsort of sorted took " << clock_selectionsort_sort << " seconds" << endl;          //prints the time it took to sort an unsorted array (seconds)

	delete[] data;                                                   //deallocate

	const size_t n2 = 70000000;                              //initialize n2 for quicksort
	size_t ARRAY_SIZE2 = n2;                                 //array size is n2
	size_t number_of_element2 = n2;
	data = new int[ARRAY_SIZE2];

	for (int i = 0; i < number_of_element2; i++) {                 //random value from 1 - ARRAY_SIZE2 into all index
		data[i] = rand() % ARRAY_SIZE2 + 1;
	}
	

	clock_t start_quicksort_unsort = clock();                        //starts clock for quicksort for unsorted

	quicksort(data, number_of_element2);                            //calls quicksort
	 
	clock_t end_quicksort_unsort = clock();                          //ends clock for quicksort for unsorted

	double clock_quicksort_unsort = double(end_quicksort_unsort - start_quicksort_unsort) / CLOCKS_PER_SEC;             //calculates the seconds it took

	cout << "Quicksort of a unsorted took " << clock_quicksort_unsort << " seconds" << endl;

	delete[] data;                                               //deallocate

	const size_t n3 = 700000000;                                 //initialize n3 for mergesort
	size_t ARRAY_SIZE3 = n3;                                     //array size is n3
	size_t number_of_element3 = n3;
	data = new int[ARRAY_SIZE3];

	for (int i = 0; i < number_of_element3; i++) {                    //random value from 1 - ARRAY_SIZE3 in each index
		data[i] = rand() % ARRAY_SIZE3 + 1;
	}

	clock_t start_mergesort_unsort = clock();                         //starts  clock for mergesort (unsort)

	mergesort(data, number_of_element3);                              //calls mergesort

	clock_t end_mergesort_unsort = clock();                           //ends clock for mergesort (unsort)

	double clock_mergesort_unsort = double(end_mergesort_unsort - start_mergesort_unsort) / CLOCKS_PER_SEC;               //converts into seconds

	cout << "Mergesort of a unsorted took " << clock_mergesort_unsort << " seconds" << endl;            //output seconds took to sort unsorted array
	
	clock_t start_mergesort_sort = clock();                                     //starts clock for mergesort (sorted)

	mergesort(data, number_of_element3);                                       //calls mergesort

	clock_t end_mergesort_sort = clock();                                       //ends clock for mergesort (sorted)

	double clock_mergesort_sort = double(end_mergesort_sort - start_mergesort_sort) / CLOCKS_PER_SEC;              //converts to seconds

	cout << "Mergesort of a sorted took " << clock_mergesort_sort << " seconds" << endl;                       //output seconds took for mergesort for sorted array

	delete[] data;                                 //deallocate

	const size_t n4 = 200000000;                      //n4 for heapsort
	size_t ARRAY_SIZE4 = n4;                          //size of array is n4
	size_t number_of_element4 = n4;
	data = new int[ARRAY_SIZE4];                       //re-allocate

	for (int i = 0; i < number_of_element4; i++) {          //random ineger from i - ARRAY_SIZE4
		data[i] = rand() % ARRAY_SIZE4 + 1;
	}

	clock_t start_heapsort_unsort = clock();                 //starts clock for heapsort (unsort)

	heapsort(data, number_of_element4);                      //calls heapsort

	clock_t end_heapsort_unsort = clock();                    //ends clock for heapsort (unsort)

	double clock_heapsort_unsort = double(end_heapsort_unsort - start_heapsort_unsort) / CLOCKS_PER_SEC;              //convert into seconds

	cout << "Heapsort of a unsort array took " << clock_heapsort_unsort << " seconds" << endl;                   //output seconds for heapsort (unsort)

	clock_t start_heapsort_sort = clock();                        //starts clock for heapsort (sorted)

	heapsort(data, number_of_element4);                           //calls heapsort

	clock_t end_heapsort_sort = clock();                         //ends clock for heapsort (sorted)

	double clock_heapsort_sort = double(end_heapsort_sort - start_heapsort_sort) / CLOCKS_PER_SEC;        //converts to seconds

	cout << "Heapsort of a sort array took " << clock_heapsort_sort << " seconds" << endl;                 //output seconds for heapsort (sorted)
	 
	delete[] data;                                     //deallocate when not in used

	return EXIT_SUCCESS;                               //end
}

void insertionsort(int data[], size_t n)                 //insertion sort
// Library facilities used: algorithm, cstdlib
{

	int i, j;
	int key;
	for (j = 1; j < n; j++)
	{
		key = data[j];
		i = j - 1;
		while (data[i] > key && i >= 0)
		{
			data[i + 1] = data[i];
			i--;
		}
		data[i + 1] = key;
	}


}

void selectionsort(int data[], size_t n)             //selection sort
// Library facilities used: algorithm, cstdlib
{

	size_t i, j, index_of_largest;
	int largest;

	if (n == 0)
		return; // No work for an empty array.

	for (i = n - 1; i > 0; --i)
	{
		largest = data[0];
		index_of_largest = 0;
		for (j = 1; j <= i; ++j)
		{
			if (data[j] > largest)
			{
				largest = data[j];
				index_of_largest = j;
			}
		}
		swap(data[i], data[index_of_largest]);
	}

}


void quicksort(int data[], size_t n)                                  //quicksort
// Library facilities used: cstdlib
{
	size_t pivot_index; // Array index for the pivot element
	size_t n1;          // Number of elements before the pivot element
	size_t n2;          // Number of elements after the pivot element

	if (n > 1)
	{
		// Partition the array, and set the pivot index.
		partition(data, n, pivot_index);

		// Compute the sizes of the subarrays.
		n1 = pivot_index;
		n2 = n - n1 - 1;

		// Recursive calls will now sort the subarrays.
		quicksort(data, n1);
		quicksort((data + pivot_index + 1), n2);
	}
}

static
void partition(int data[], size_t n, size_t& pivot_index)
// Library facilities used: itemtool.h, stdlib.h
//
// NOTES FROM THE IMPLEMENTOR:
// How the partition works on small arrays:
//
// Notice that n=0 is not permitted by the precondition.
//
// If n=1, then too_big_index is initialized as 1, and too_small_index is
// initialized as 0. Therefore, the body of the loop is never executed,
// and after the loop pivot_index is set to zero.
//
// If n=2, then both too_big_index and too_small_index are initialized as 1.
// The loop is entered, and there are two cases to consider:
// -- if data[1] <= pivot, then too_big_index increases to 2, and
//    too_small_index stays at 1. The if-statement at the bottom of the loop
//    is then skipped, and after the loop we copy data[1] down to data[0],
//    and copy the pivot into data[0]. Thus, the smaller element is in
//    data[0], and the larger element (the pivot) is in data[1].
// -- if data[1] > pivot, then too_big_index stays at 1, and too_small_index
//    decreases to 0. The if-statement at the bottom of the loop is then
//    skipped, and after the loop we end up copying the pivot onto data[0]
//    (leaving data[1] alone). Thus, the smaller element (the pivot) remains
//    at data[0], leaving the larger element at data[1].
{
	int pivot = data[0];
	size_t too_big_index = 1;     // Index of first item after pivot
	size_t too_small_index = n - 1; // Index of last item

	// Partition the array, using pivot as the pivot element
	while (too_big_index <= too_small_index)
	{
		while ((too_big_index < n) && (data[too_big_index] <= pivot))
			too_big_index++;
		while (data[too_small_index] > pivot)
			too_small_index--;
		if (too_big_index < too_small_index)
			swap(data[too_small_index], data[too_big_index]);
	}

	// Move the pivot element to its correct position
	pivot_index = too_small_index;
	data[0] = data[pivot_index];
	data[pivot_index] = pivot;
}

void mergesort(int data[], size_t n)
// Precondition: data is an array with at least n components.
// Postcondition: The elements of data have been rearranged so
// that data[0] <= data[1] <= ... <= data[n-1].
// NOTE: If there is insufficient dynamic memory, thenbad_alloc is thrown.
// Library facilities used: cstdlib
{
	size_t n1; // Size of the first subarray
	size_t n2; // Size of the second subarray

	if (n > 1)
	{
		// Compute sizes of the subarrays.
		n1 = n / 2;
		n2 = n - n1;

		mergesort(data, n1);         // Sort from data[0] through data[n1-1]
		mergesort((data + n1), n2);  // Sort from data[n1] to the end

		// Merge the two sorted halves.
		merge(data, n1, n2);
	}
}

void merge(int data[], size_t n1, size_t n2)                                           //mergesort
// Precondition: data is an array (or subarray) with at least n1 + n2 elements.
// The first n1 elements (from data[0] to data[n1 - 1]) are sorted from
// smallest to largest, and the last n2 (from data[n1] to data[n1 + n2 - 1])
// also are sorted from smallest to largest.
// Postcondition: The first n1 + n2 elements of data have been rearranged to be
// sorted from smallest to largest.
// NOTE: If there is insufficient dynamic memory, then bad_alloc is thrown.
// Library facilities used: cstdlib
{
	int* temp;          // Points to dynamic array to hold the sorted elements
	size_t copied = 0; // Number of elements copied from data to temp
	size_t copied1 = 0; // Number copied from the first half of data
	size_t copied2 = 0; // Number copied from the second half of data
	size_t i;           // Array index to copy from temp back into data

	// Allocate memory for the temporary dynamic array.
	temp = new int[n1 + n2];

	// Merge elements, copying from two halves of data to the temporary array.
	while ((copied1 < n1) && (copied2 < n2))
	{
		if (data[copied1] < (data + n1)[copied2])
			temp[copied++] = data[copied1++];        // Copy from first half
		else
			temp[copied++] = (data + n1)[copied2++]; // Copy from second half
	}

	// Copy any remaining entries in the left and right subarrays.
	while (copied1 < n1)
		temp[copied++] = data[copied1++];
	while (copied2 < n2)
		temp[copied++] = (data + n1)[copied2++];

	// Copy from temp back to the data array, and release temp's memory.
	for (i = 0; i < n1 + n2; i++)
		data[i] = temp[i];
	delete[] temp;
}

void heapsort(int data[], size_t n)                                        //heapsort
// Library facilities used: algorithm, cstdlib
{
	size_t unsorted;

	make_heap(data, n);

	unsorted = n;

	while (unsorted > 1)
	{
		--unsorted;
		swap(data[0], data[unsorted]);
		reheapify_down(data, unsorted);
	}
}

size_t parent(size_t k)
// Library facilities used: cstdlib
{
	return (k - 1) / 2;
}

size_t left_child(size_t k)
// Library facilities used: cstdlib
{
	return 2 * k + 1;
}

size_t right_child(size_t k)
// Library facilities used: cstdlib
{
	return 2 * k + 2;
}

void make_heap(int data[], size_t n)
// Library facilities used: itemtool.h (from page 277), cstdlib
// 
{
	size_t i;  // Index of next element to be added to heap
	size_t k;  // Index of new element as it is being pushed upward through the heap

	for (i = 1; i < n; ++i)
	{
		k = i;
		while ((k > 0) && (data[k] > data[parent(k)]))
		{
			swap(data[k], data[parent(k)]);
			k = parent(k);
		}
	}
}

void reheapify_down(int data[], size_t n)
// Library facilities used: itemtool.h (from page 277), cstdlib
{
	size_t current;          // Index of the node that's moving down
	size_t big_child_index;  // Index of the larger child of the node that's moving down
	bool heap_ok = false;    // Will change to true when the heap becomes correct

	current = 0;

	// Note: The loop keeps going while the heap is not okay, and while the current node has
	// at least a left child. The test to see whether the current node has a left child is
	// left_child(current) < n.
	while ((!heap_ok) && (left_child(current) < n))
	{
		// Compute the index of the larger child:
		if (right_child(current) >= n)
			// There is no right child, so left child must be largest
			big_child_index = left_child(current);
		else if (data[left_child(current)] > data[right_child(current)])
			// The left child is the bigger of the two children
			big_child_index = left_child(current);
		else
			// The right child is the bigger of the two children
			big_child_index = right_child(current);

		// Check whether the larger child is bigger than the current node. If so, then swap
		// the current node with its bigger child and continue; otherwise we are finished.
		if (data[current] < data[big_child_index])
		{
			swap(data[current], data[big_child_index]);
			current = big_child_index;
		}
		else
			heap_ok = true;
	}
}
