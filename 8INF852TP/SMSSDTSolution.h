
#pragma once
#ifndef __SMSSDTSOLUTION_H
#define __SMSSDTSOLUTION_H
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "SMSSDTProblem.h"

using namespace std;

class SMSSDTSolution
{
protected:
	double FctObj;

public:

	vector<int > Solution;
	vector<int > CT;
	vector<int > TT;
	vector<int > ST;


	/**
	 * SMSSDTSolution()
	 * Constructeur de la classe
	**/
	SMSSDTSolution(int N) {
		Solution.resize(N);
		TT.resize(N);
		CT.resize(N);
		ST.resize(N);
		FctObj = -1;

	};

	/**
	 * SMSSDTSolution(int N)
	 * Constructeur de la classe
	 * Assigne les valeurs al�atoires initiales � la solution
	**/
	SMSSDTSolution(int N, bool test);


	/**
	 * SMSSDTSolution(SMSSDTProblem* LeProb, SMSSDTSolution* Sol);
	 * Constructeur de la classe
	 * g�n�rer la meilleure solution voisinnage � Sol en utilisant Descente
	**/
	SMSSDTSolution(SMSSDTProblem* LeProb, SMSSDTSolution p1, SMSSDTSolution p2);
	SMSSDTSolution(SMSSDTProblem* LeProb, SMSSDTSolution& Sol);
	SMSSDTSolution(SMSSDTProblem* LeProb);
	/* Destructeur*/
	~SMSSDTSolution();


	/**
	 * Save
	 * Fonction permettant d'�crire la solution dans un flux
	 * @param Stream : Flux dans lequel on doit �crire
	**/
	ostream& Save(ostream& Stream);


	/**  definition des accesseurs **/

	/* Retourne la valeur de l'objectif */
	double getObj() {
		return FctObj;
	}

	void  inverse(int N, int a, int b);
	void opt(int N, int a, int b);
	void setObj(double obj) {
		FctObj = obj;
	}
};




#endif