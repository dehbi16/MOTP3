#ifndef __SMSSDTPROBLEM_CPP
#define __SMSSDTPROBLEM_CPP


#include "SMSSDTSolution.h"
#include "Tools.h"

#include <cmath>
#include <conio.h>
using namespace std;

/**
 * SMSSDTSolution()
 * Constructeur de la classe
 * Assigne les valeurs aléatoires initiales à la solution
**/

SMSSDTSolution::SMSSDTSolution(int N, bool test) {
	// scheduling vector
	Solution.resize(N);
	TT.resize(N);
	CT.resize(N);
	ST.resize(N);
	// initialisation of possible values
	vector < int > possibles(N);
	for (int i = 0; i < N; i++)
		possibles[i] = i;
	// random initialization
	int rInd;              // random index
	for (int i = 0; i < N; i++)
	{
		rInd = (int)(((double)rand() / ((double)RAND_MAX + 1.0)) * (N - i));//rng.uniform (N - i);
		Solution[i] = possibles[rInd];
		possibles[rInd] = possibles[N - i - 1];
	}
}

SMSSDTSolution::SMSSDTSolution(SMSSDTProblem* LeProb, SMSSDTSolution& Sol) {
	this->setObj(-1);
	int N = LeProb->getN();
	Solution.resize(N);
	TT.resize(N);
	CT.resize(N);
	ST.resize(N);
	vector<int> s = Sol.Solution;
	Tools::swapMove(Sol.Solution, Solution, N);

}


SMSSDTSolution::SMSSDTSolution(SMSSDTProblem* LeProb, SMSSDTSolution p1, SMSSDTSolution p2) {
	this->setObj(-1);
	int N = LeProb->getN();
	Solution.resize(N);
	TT.resize(N);
	CT.resize(N);
	ST.resize(N);
	vector<int>  ec(N);
	vector<int> P1 = p1.Solution;
	vector<int> P2 = p1.Solution;

	vector <vector < int > > s = LeProb->getS();
	int rInd = (int)(((double)rand() / ((double)RAND_MAX + 1.0)) * (N - 1));
	int index1 = 0;
	int a, b;

	if (rInd < N / 2) {
		Solution[0] = P1[0];
		ec[P1[0]] = 1;
	}
	else { Solution[0] = P2[0]; ec[P2[0]] = 1; }

	for (int i = 1; i < N; i++) {
		a = P1[i]; b = P2[i];
		if (ec[P1[i]] == 1) {
			while (ec[index1] == 1) index1++;
			a = index1;
			index1 = 0;
		}
		if (ec[P2[i]] == 1) {
			while (ec[index1] == 1) index1++;
			b = index1;
			index1 = 0;
		}

		if (s[i - 1][a] < s[i - 1][b]) {
			Solution[i] = a;
			ec[a] = 1;
		}
		else {
			Solution[i] = b;
			ec[b] = 1;
		}
	}


}
// CONSTRUCTEUR REGLE EDD
SMSSDTSolution::SMSSDTSolution(SMSSDTProblem* LeProb) {
	this->setObj(-1);
	int N = LeProb->getN();
	std::vector < int > d(N), p(N);
	d = LeProb->getD();
	Solution.resize(N);
	TT.resize(N);
	CT.resize(N);
	ST.resize(N);
	for (int k = 0; k < N; k++) Solution[k] = k;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (d[Solution[i]] > d[Solution[j]]) {
				int x = Solution[i];
				Solution[i] = Solution[j];
				Solution[j] = x;
			}
		}
	}
}

/**
 * Save
 * Fonction permettant d'écrire la solution dans un flux
 * @param Stream : Flux dans lequel on doit écrire
**/
ostream& SMSSDTSolution::Save(ostream& Stream) {

	for (int i = 0; i < (int)Solution.size(); i++) {
		Stream << Solution[i] << " ";

	}

	return Stream;
}


SMSSDTSolution::~SMSSDTSolution() {
	Solution.clear();
}




void SMSSDTSolution::opt(int N, int a, int b) {



	int x = Solution[a];
	Solution[a] = Solution[b];
	Solution[b] = x;
	//inverse(N, a, b);

}

#endif  //__SMSSDTPROBLEM_CPP