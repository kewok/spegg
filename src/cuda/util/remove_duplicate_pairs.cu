#include <util/remove_duplicate_pairs.h>

#include <thrust/copy.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/host_vector.h>

// This code is inspired by: https://stackoverflow.com/questions/5521091/thrust-removing-duplicates-in-key-value-arrays

void remove_duplicate_pairs(thrust::device_vector<int> &vectorA, 
			    thrust::device_vector<int> &vectorB )
	{
	/*
	************************
	//
	//
	// Helper function for identifying and removing duplicates in pairs made up of two vectors; e.g.
	// vectorA = [0,0,0, 1,1,1, 2,2,2, 0,0]
	// vectorB = [3,4,5, 6,6,7, 8,9,9, 6,3]
	// 
	// Remove pairs so that:
	//
	// vectorA = [0,0,0,0, 1,1, 2,2]
	// vectorB = [3,4,5,6, 6,7, 8,9]
	//
	// Note that vectorA and vectorB must be of the same size.
	//
	************************
	*/

	if (vectorA.size() != vectorB.size())
		{
		std::cout << "vectors in remove_duplicate_pairs have to be the same size!" << std::endl;
		return;
		}

	int total_vector_sizes_before_removing_duplicates = vectorA.size(); 
	// Check for duplicates 
	thrust::device_vector<int> copy_vectorA( total_vector_sizes_before_removing_duplicates );
	thrust::device_vector<int> copy_vectorB( total_vector_sizes_before_removing_duplicates );
	
	thrust::copy(vectorA.begin(), vectorA.end(), copy_vectorA.begin());
	thrust::copy(vectorB.begin(), vectorB.end(), copy_vectorB.begin()); 

	// Sort by vectorA and then by vectorB (e.g.,vectorA = [0,0,0, 0,1,1, 2,2,2, 0,0], vectorB = [3,4,5, 6,6,7, 8,9,9, 6,3] becomes vectorA = [0,0,0,0,0,0,1,1,2,2,2], vectorB = [3,3,4,5,6,6,6,7,8,9,9] )

	thrust::sort(thrust::make_zip_iterator(thrust::make_tuple(copy_vectorA.begin(), copy_vectorB.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(copy_vectorA.end(), copy_vectorB.end())));

	// Remove duplicates, based on contiguity in vectorB

	typedef thrust::device_vector< int >                IntVector;
	typedef IntVector::iterator                         IntIterator;
	typedef thrust::tuple< IntIterator, IntIterator >   IntIteratorTuple;
	typedef thrust::zip_iterator< IntIteratorTuple >    ZipIterator;

	ZipIterator newEnd = thrust::unique( thrust::make_zip_iterator( thrust::make_tuple( copy_vectorA.begin(), copy_vectorB.begin() ) ), thrust::make_zip_iterator( thrust::make_tuple( copy_vectorA.end(), copy_vectorB.end() ) ) );

	IntIteratorTuple endTuple = newEnd.get_iterator_tuple();

// by magic, this makes sure you have the correct length of copy_vectorA, which should be copy_vectorA.size() - duplicates:

	copy_vectorA.erase( thrust::get<0>( endTuple ),copy_vectorA.end() );
	copy_vectorB.erase( thrust::get<1>( endTuple ), copy_vectorB.end() );

	// Correct the values for vectorA and vectorB to remove the duplicates

	vectorA.resize( copy_vectorA.size() );
	vectorB.resize( copy_vectorA.size() );

	thrust::fill(vectorA.begin(), vectorA.end(), 0);
	thrust::fill(vectorB.begin(),vectorB.end(), 0);

	thrust::copy(copy_vectorA.begin(), copy_vectorA.end(), vectorA.begin());
	thrust::copy(copy_vectorB.begin(), copy_vectorB.end(), vectorB.begin());
	}

