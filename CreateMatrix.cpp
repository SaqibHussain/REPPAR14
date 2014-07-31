/***************************************************
File:			CreateMatrix.cpp
Description:		Random matrix generator
Author:			Saqib Hussain 
Organisation:		School of Engineering, Cranfield University
Email:			s.m.hussain@cranfield.ac.uk
Copyright:		Copyright Saqib Hussain 2014
****************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
	int N = atoi(argv[1]);
	int M = atoi(argv[2]);
	double** mat = new double*[N];
	for(int i = 0; i < N; ++i)
		mat[i] = new double[M];

	for (int i=0; i<N; i++) {
		for (int j=0; j<M; j++)
			mat[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	}

	ofstream ofs;
	ofs.open("data.txt");

	ofs << N << " " << M << "\n";
	for (int i=0; i<N; i++) {
		for (int j=0; j<M; j++)
			ofs << mat[i][j] << " ";
		ofs << "\n";
	}
}
