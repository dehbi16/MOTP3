
#include <iostream>
#include <tuple>
#include "Tools.h"
#include <algorithm>  
using namespace std;

void mutation(vector<SMSSDTSolution*>* SolutionsC, int N) {
	double a;
	int index1; int index2;
	for (int i = 0; i < SolutionsC->size(); i++) {
		a = ((double)rand()) / RAND_MAX;
		if (a > 0.66) {

			index1 = rand() % N;
			do {
				index2 = rand() % N;
			} while (index1 == index2);

			(*SolutionsC)[i]->Solution[index1] = (*SolutionsC)[i]->Solution[index1] + (*SolutionsC)[i]->Solution[index2];
			(*SolutionsC)[i]->Solution[index2] = (*SolutionsC)[i]->Solution[index1] - (*SolutionsC)[i]->Solution[index2];
			(*SolutionsC)[i]->Solution[index1] = (*SolutionsC)[i]->Solution[index1] - (*SolutionsC)[i]->Solution[index2];
			/*
			Tools::scrambleMove((*SolutionsC)[i]->Solution, (*SolutionsC)[i]->Solution, N, 4);
			*/
		}
	}

}

void crossOver(vector<SMSSDTSolution*>* SolutionsC, vector<int> P1, vector<int> P2) {
	int const N = P1.size();
	SMSSDTSolution* solutionC1 = new SMSSDTSolution(N);
	SMSSDTSolution* solutionC2 = new SMSSDTSolution(N);

	vector<int> E1;
	vector<int> E2;
	E1.resize(N);
	E2.resize(N);

	vector<int>  e1check(N);
	vector<int>  e2check(N);

	// choisir deux indices

	int i = (int)rand() % N;
	int j;
	do {
		j = (int)rand() % N;
	} while (i == j);



	// i doit tjrs moins que j
	if (i > j) {
		i = i + j;
		j = i - j;
		i = i - j;
	}

	//int i = N / 4;
	//int j = (3 * N) / 4; */


	int index1 = 0;
	int index2 = 0;
	int k;

	for (k = 0; k < i; k++) {
		if (e1check[P1[k]] == 1) {
			while (e1check[index1] == 1) index1++;
			E1[k] = index1;
			e1check[index1] = 1;
		}
		else {
			E1[k] = P1[k];
			e1check[P1[k]] = 1;
		}
		if (e2check[P2[k]] == 1) {
			while (e2check[index2] == 1) index2++;
			E2[k] = index2;
			e2check[index2] = 1;
		}
		else {
			E2[k] = P2[k];
			e2check[P2[k]] = 1;
		}
	}

	for (k = i; k < j; k++) {
		if (e1check[P2[k]] == 1) {
			while (e1check[index1] == 1) index1++;
			E1[k] = index1;
			e1check[index1] = 1;
		}
		else {
			E1[k] = P2[k];
			e1check[P2[k]] = 1;
		}
		if (e2check[P1[k]] == 1) {
			while (e2check[index2] == 1) index2++;
			E2[k] = index2;
			e2check[index2] = 1;
		}
		else {
			E2[k] = P1[k];
			e2check[P1[k]] = 1;
		}
	}
	for (k = j; k < N; k++) {
		if (e1check[P1[k]] == 1) {
			while (e1check[index1] == 1) index1++;
			E1[k] = index1;
			e1check[index1] = 1;
		}
		else {
			E1[k] = P1[k];
			e1check[P1[k]] = 1;
		}
		if (e2check[P2[k]] == 1) {
			while (e2check[index2] == 1) index2++;
			E2[k] = index2;
			e2check[index2] = 1;
		}
		else {
			E2[k] = P2[k];
			e2check[P2[k]] = 1;
		}
	}

	solutionC1->Solution = E1;
	solutionC2->Solution = E2;

	SolutionsC->push_back(solutionC1);
	SolutionsC->push_back(solutionC2);
}

void crossOver2(SMSSDTProblem* LeProb, vector<SMSSDTSolution*>* SolutionsC, vector<int> P1, vector<int> P2) {
	int const N = P1.size();
	SMSSDTSolution* solutionC1 = new SMSSDTSolution(N);
	SMSSDTSolution* solutionC2 = new SMSSDTSolution(N);

	vector<int> E1;
	vector<int> E2;
	E1.resize(N);
	E2.resize(N);

	vector<int>  e1check(N);
	vector<int>  e2check(N);

	// choisir deux indices

	int i = (int)rand() % N;
	int j;
	do {
		j = (int)rand() % N;
	} while (i == j);



	// i doit tjrs moins que j
	if (i > j) {
		i = i + j;
		j = i - j;
		i = i - j;
	}

	//int i = N / 4;
	//int j = (3 * N) / 4; */


	int index1 = 0;
	int index2 = 0;
	int k;

	for (k = 0; k < i; k++) {
		if (e1check[P1[k]] == 1) {
			while (e1check[index1] == 1) index1++;
			E1[k] = index1;
			e1check[index1] = 1;
			index1 = 0;
		}
		else {
			E1[k] = P1[k];
			e1check[P1[k]] = 1;
		}
		if (e2check[P2[k]] == 1) {
			while (e2check[index2] == 1) index2++;
			E2[k] = index2;
			e2check[index2] = 1;
			index2 = 0;
		}
		else {
			E2[k] = P2[k];
			e2check[P2[k]] = 1;
		}
	}
	int a = 0; int b = 0;
	vector <vector < int > > s = LeProb->getS();


	for (k = i; k < j - 1; k++) {
		a += s[P1[i]][P1[i + 1]];
		b += s[P2[i]][P2[i + 1]];
	}
	if (a < b) {
		for (k = i; k < j; k++) {
			if (e1check[P1[k]] == 1) {
				while (e1check[index1] == 1) index1++;
				E1[k] = index1;
				e1check[index1] = 1;
				index1 = 0;
			}
			else {
				E1[k] = P1[k];
				e1check[P1[k]] = 1;
			}


			if (e2check[P1[k]] == 1) {
				while (e2check[index2] == 1) index2++;
				E2[k] = index2;
				e2check[index2] = 1;
				index2 = 0;
			}
			else {
				E2[k] = P1[k];
				e2check[P1[k]] = 1;
			}
		}
	}
	else {


		for (k = i; k < j; k++) {
			if (e1check[P2[k]] == 1) {
				while (e1check[index1] == 1) index1++;
				E1[k] = index1;
				e1check[index1] = 1;
				index1 = 0;
			}
			else {
				E1[k] = P2[k];
				e1check[P2[k]] = 1;
			}


			if (e2check[P2[k]] == 1) {
				while (e2check[index2] == 1) index2++;
				E2[k] = index2;
				e2check[index2] = 1;
				index2 = 0;
			}
			else {
				E2[k] = P2[k];
				e2check[P2[k]] = 1;
			}
		}
	}



	for (k = j; k < N; k++) {
		if (e1check[P2[k]] == 1) {
			while (e1check[index1] == 1) index1++;
			E1[k] = index1;
			e1check[index1] = 1;
			index1 = 0;
		}
		else {
			E1[k] = P2[k];
			e1check[P2[k]] = 1;
		}
		if (e2check[P1[k]] == 1) {
			while (e2check[index2] == 1) index2++;
			E2[k] = index2;
			e2check[index2] = 1;
			index1 = 0;
		}
		else {
			E2[k] = P1[k];
			e2check[P1[k]] = 1;
		}
	}

	solutionC1->Solution = E1;
	solutionC2->Solution = E2;

	SolutionsC->push_back(solutionC1);
	SolutionsC->push_back(solutionC2);
}

void crossOver3(vector<SMSSDTSolution*>* SolutionsC, vector<int> P1, vector<int> P2) {
	int const N = P1.size();
	SMSSDTSolution* solutionC1 = new SMSSDTSolution(N);
	SMSSDTSolution* solutionC2 = new SMSSDTSolution(N);

	vector<int> E1;
	vector<int> E2;
	E1.resize(N);
	E2.resize(N);

	vector<int>  e1check(N);
	vector<int>  e2check(N);

	// choisir deux indices

	int i = (int)rand() % N;
	int j;
	do {
		j = (int)rand() % N;
	} while (i == j);



	// i doit tjrs moins que j
	if (i > j) {
		i = i + j;
		j = i - j;
		i = i - j;
	}
	/*
	int i = N / 4;
	int j = (3 * N) / 4;
	*/
	int k;
	for (k = 0; k < i; k++) {
		if (e1check[P1[k]] == 1) {
			E1[k] = -1;
		}
		else {
			E1[k] = P1[k];
			e1check[P1[k]] = 1;
		}
		if (e2check[P2[k]] == 1) {
			E2[k] = -1;
		}
		else {
			E2[k] = P2[k];
			e2check[P2[k]] = 1;
		}
	}

	for (k = i; k < j; k++) {
		if (e1check[P2[k]] == 1) {
			E1[k] = -1;
		}
		else {
			E1[k] = P2[k];
			e1check[P2[k]] = 1;
		}
		if (e2check[P1[k]] == 1) {
			E2[k] = -1;
		}
		else {
			E2[k] = P1[k];
			e2check[P1[k]] = 1;
		}
	}

	for (k = j; k < N; k++) {
		if (e1check[P1[k]] == 1) {
			E1[k] = -1;
		}
		else {
			E1[k] = P1[k];
			e1check[P1[k]] = 1;
		}
		if (e2check[P2[k]] == 1) {
			E2[k] = -1;
		}
		else {
			E2[k] = P2[k];
			e2check[P2[k]] = 1;
		}
	}
	int index1 = 0;
	int index2 = 0;
	for (k = 0; k < N; k++) {
		if (E1[k] == -1) {
			while (e1check[index1] != 0) {
				index1++;
			}
			E1[k] = index1;
			e1check[index1] = 1;
		}
		if (E2[k] == -1) {
			while (e2check[index2] != 0) {
				index2++;
			}
			E2[k] = index2;
			e2check[index2] = 1;
		}
	}


	solutionC1->Solution = E1;
	solutionC2->Solution = E2;

	SolutionsC->push_back(solutionC1);
	SolutionsC->push_back(solutionC2);
}

void voi(SMSSDTSolution& Sol, SMSSDTProblem* LeProb) { // fonction qui parcours le voisinage 
	double	dTheBestFitness = 100000;
	SMSSDTSolution	Smeilleur(LeProb->getN());

	for (int i = 0; i < LeProb->getN(); i++) {

		for (int j = i + 1; j < LeProb->getN(); j++) {

			Sol.opt(LeProb->getN(), i, j);


			Tools::Evaluer(LeProb, Sol);

			if (Sol.getObj() < dTheBestFitness) // Si améliore meilleure solution, la garder
			{
				Smeilleur = Sol;
				dTheBestFitness = Smeilleur.getObj();

			}
			Sol.opt(LeProb->getN(), i, j);



		}
	}
	Sol = Smeilleur;
}

void voi2(SMSSDTSolution& Sol, SMSSDTProblem* LeProb) { // fonction qui parcours le voisinage 
	double	dTheBestFitness = 100000;
	SMSSDTSolution	Smeilleur(LeProb->getN());
	SMSSDTSolution	voisins = (LeProb->getN());

	for (int i = 0; i < LeProb->getN() - 1; i++) {

		for (int j = i + 2; j < LeProb->getN() - 1; j++) {
			voisins = Sol;

			voisins.Solution[j] = Sol.Solution[i];
			voisins.Solution[j + 1] = Sol.Solution[i + 1];

			for (int k = i + 2; k < j + 2; k++) {

				voisins.Solution[k - 2] = Sol.Solution[k];


			}





			Tools::Evaluer(LeProb, voisins);

			if (voisins.getObj() < dTheBestFitness) // Si améliore meilleure solution, la garder
			{
				Smeilleur = voisins;
				dTheBestFitness = Smeilleur.getObj();

			}




		}
	}
	Sol = Smeilleur;
}

void voi3(SMSSDTSolution& Sol, SMSSDTProblem* LeProb) { // fonction qui parcours le voisinage 
	double	dTheBestFitness = 100000;
	SMSSDTSolution	Smeilleur(LeProb->getN());
	SMSSDTSolution	voisins = (LeProb->getN());

	for (int i = 0; i < LeProb->getN() - 2; i++) {

		for (int j = i + 3; j < LeProb->getN() - 2; j++) {
			voisins = Sol;

			voisins.Solution[j] = Sol.Solution[i];
			voisins.Solution[j + 1] = Sol.Solution[i + 1];

			voisins.Solution[j + 2] = Sol.Solution[i + 2];
			for (int k = i + 3; k < j + 3; k++) {

				voisins.Solution[k - 2] = Sol.Solution[k];


			}





			Tools::Evaluer(LeProb, voisins);

			if (voisins.getObj() < dTheBestFitness) // Si améliore meilleure solution, la garder
			{
				Smeilleur = voisins;
				dTheBestFitness = Smeilleur.getObj();

			}




		}
	}
	Sol = Smeilleur;
}

void desente(SMSSDTProblem* LeProb, SMSSDTSolution* pSolution, SMSSDTSolution& Smeilleur) {

	int N = LeProb->getN();
	//sl.resize(N);
	double	dTheBestFitness = 100000;	//Fitness de la meilleure solution
	int e = 0;
	Tools::Evaluer(LeProb, *pSolution);	//Évaluer la solution
	dTheBestFitness = pSolution->getObj();
	while (e == 0) {
		SMSSDTSolution	sol(LeProb->getN());
		sol = *pSolution;
		voi(sol, LeProb);
		if (sol.getObj() < dTheBestFitness) // Si la meilleur sol du voisinage  améliore meilleure solution, la garder
		{
			*pSolution = sol;
			dTheBestFitness = pSolution->getObj();
			//sl.push_back(pSolution->getObj());
		}
		else { e = 1; }
	}
	Smeilleur = *pSolution;

}

void shaking(SMSSDTProblem* LeProb, vector<int> const& solution, vector<int>& voisins, int mode) {
	switch (mode) {
	case 0:
		Tools::swapMove(solution, voisins, LeProb->getN());
		break;
	case 1:
		Tools::insertionMove(solution, voisins, LeProb->getN());
		break;
	case 2:

		Tools::EDDMove(solution, voisins, LeProb->getD(), LeProb->getN(), rand() % (LeProb->getN()));
		break;
	case 4:

		Tools::scrambleMove(solution, voisins, LeProb->getN(), rand() % (LeProb->getN()));
		break;
	case 5:
		Tools::inversionMove(solution, voisins, LeProb->getN(), rand() % (LeProb->getN()));

		break;
	}
}

float p(float T, double x, double xp) {
	return exp((x - xp) / T);
}

float g(float T) {
	return 0.75 * T;
}

void recuitSimule(SMSSDTProblem* LeProb, SMSSDTSolution* solution, SMSSDTSolution& resultat) {
	Tools::Evaluer(LeProb, *solution);
	SMSSDTSolution* pSolution;
	float T = 20;
	float r;
	for (int i = 0; i < 200; i++) {
		pSolution = new SMSSDTSolution(LeProb, *solution);
		Tools::Evaluer(LeProb, *pSolution);	//Évaluer la solution
		r = (float)rand() / (float)RAND_MAX;
		if (r < p(T, solution->getObj(), pSolution->getObj())) {
			solution = pSolution;
		}
		T = g(T);
	}
	resultat = *solution;
}

void vns(SMSSDTSolution* pSolution, SMSSDTSolution& Smeilleur, SMSSDTProblem* LeProb) {
	SMSSDTSolution	Svoisin = NULL;
	SMSSDTSolution	Svoisin1 = NULL;
	Svoisin = *pSolution;
	Svoisin1 = *pSolution;
	int amelioration = 0;
	int m = 0;
	while (amelioration < 50) {
		shaking(LeProb, pSolution->Solution, Svoisin.Solution, m);
		desente(LeProb, &Svoisin, Svoisin1);

		if (Svoisin1.getObj() < pSolution->getObj()) {
			*pSolution = Svoisin1;
			m = -1;

			amelioration = 0;
		}
		else {
			amelioration++;
		}
		if (m < 5) {
			m++;
		}
		else { m = 0; }
		if (pSolution->getObj() == 0) {
			break;
		}
	}
	amelioration = 0;
	Smeilleur = *pSolution;
}

void tabou(SMSSDTSolution* pSolution, SMSSDTSolution& Smeilleur, SMSSDTProblem* LeProb) {
	vector<vector <int>> Ta;
	int N = 6;
	int index = 0;
	Ta.resize(N);
	//Sauvegarde de la meilleure solution
	SMSSDTSolution	Svoisin = NULL;
	SMSSDTSolution	Svoisin1 = NULL;
	for (int i = 0; i < N; i++) Ta[i].resize(LeProb->getN());

	// INITIALISATION DE LA SOLUTION
	Svoisin = SMSSDTSolution(LeProb->getN(), true);
	Svoisin = *pSolution;
	Smeilleur = Svoisin;

	for (int i = 0; i < 200; i++) {
		// TROUVER LA SOLUTION MINIMISE LA FONCTION DANS Nt(X)
		// vérifier si Svoisin dans la liste Tabou
		do {
			Svoisin1 = SMSSDTSolution(LeProb, Svoisin);
		} while (Tools::contains(Ta, Svoisin1.Solution));

		Tools::Evaluer(LeProb, Svoisin1);
		for (int j = 0; j < 7; j++) {
			do {
				pSolution = new SMSSDTSolution(LeProb, Svoisin);
			} while (Tools::contains(Ta, pSolution->Solution));

			Tools::Evaluer(LeProb, *pSolution);
			if (pSolution->getObj() < Svoisin1.getObj()) {
				// vérifier si psolution dans la liste Tabou
				Svoisin1 = *pSolution;
			}
		}

		if (Svoisin1.getObj() < Smeilleur.getObj()) {
			Smeilleur = Svoisin1;
		}
		Svoisin = Svoisin1;

		Ta[index % N] = Svoisin.Solution;

		index++;
		/*
		for (int k = 0; k < N; k++) {
			for (int l = 0; l < LeProb->getN(); l++) {
				cout << Ta[k][l] << " ";
			}
			cout << endl;
		}
		cout<< endl;
		*/

	}

}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	clock_t	Start, End;	//Déclaration de variable afin de calculer le temps écoulé
	double Elapsed = 0;	//Variable servant à calculer le temps écoulé (Différence entre End et Start
	double	dTheBestFitness = 100000;	//Fitness de la meilleure solution
	SMSSDTProblem* LeProb;	//Déclaration d'un problème	
	LeProb = new SMSSDTProblem(argv[2]);	//Lecture du deuxi;eme paramètre à partir de la console
	//LeProb->printOn(cout);	// Imprimer le Problème
	SMSSDTSolution* pSolution = NULL;	//Solution intermédiaire


	// argv[1] exécutions de la génération aléatoire
	for (int j = 0; j < atoi(argv[1]); j++)
	{
		Start = clock();	//Démarrer l'horloge	
		SMSSDTSolution	Smeilleur(LeProb->getN());	//Sauvegarde de la meilleure solution
		
		// initialiser les paramètres
		int nbiter = 300;
		int m = 20;
		int tau0 = 1;
		int N = LeProb->getN();
		float alpha = 6;
		float beta = 3;
		int Q = 100;
		vector<SMSSDTSolution*> ants;
		ants.resize(m);


		int index;

		// Initialiser le pheromone
		vector < vector < double> > tau;
		tau.resize(N);
		for (int i = 0; i < N; i++) tau[i].resize(N);
		for (int k = 0; k < N; k++) {
			for (int l = 0; l < N; l++) {
				tau[k][l] = tau0;
			}

		}

		for (int t = 0; t < nbiter; t++) {
			// choisir le premier job pour chaque fourmis 
			vector<double> P;
			P.resize(N);
			for (int l = 0; l < m; l++) {
				pSolution = new SMSSDTSolution(N);
				index = (int)rand() % N;
				pSolution->Solution[0] = index;
				ants[l] = pSolution;
				
			}
			
			for (int k = 1; k < N; k++) {
				for (int l = 0; l < m; l++) {
					for (int i = 0; i < N; i++) P[i] = -1.0;
					for (int i = 0; i < k; i++) {
						P[ants[l]->Solution[i]] = 0.0;
					}
					
					float deno = 0;
					for (int i = 0; i < N; i++) {
						if (P[i] != 0) {
							P[i] = pow(tau[ants[l]->Solution[k - 1]][i],alpha) * pow(1.0/(LeProb->getS()[ants[l]->Solution[k - 1]][i] + LeProb->getP()[i]),beta);
							deno += P[i];
						}
					}
					for (int i = 0; i < N; i++) P[i] = P[i]/deno;
					
					//calculer la meilleur probabilité
					float max = P[0];
					index = 0;
					for (int i = 1; i < N; i++) {
						if (max < P[i]) {
							max = P[i];
							index = i;
						}
					}

					ants[l]->Solution[k] = index;
					
					
					
				}
			}

			for (int l = 0; l < m; l++) {
				desente(LeProb, ants[l], *ants[l]);
			}

			// dépose une quantité de phéromone 
			
			int min;
			for (int l = 0; l < m; l++) {
				Tools::Evaluer(LeProb, *ants[l]);
				if (l == 0) {
					min = ants[l]->getObj();
					index = l;
				}
				else {
					if (min > ants[l]->getObj()) {
						min = ants[l]->getObj();
						index = l;
					}
				}
			}
			for (int i = 1; i < N; i++) {
				tau[ants[index]->Solution[i - 1]][ants[index]->Solution[i]] += Q / ants[index]->getObj();
			}
			
			
			// Evaporation 

		}
		
		for (int i = 0; i < N; i++) {
			for (int k = 0; k < N; k++) {
				cout << tau[i][k] << " ";
			}
			cout << endl;
		}

		for (int l = 0; l < m; l++) {
			for (int i = 0; i < N; i++) {
				cout << ants[l]->Solution[i] << " ";
			}
			cout << "\t" << ants[l]->getObj() << endl;
		}
		cout << endl;


		Smeilleur = *ants[index];
		
		End = clock(); // Arrêter le clock
		Elapsed = (double(End - Start)) / CLOCKS_PER_SEC;	//Calculer le temps écoulé
		Tools::WriteReportLog(Elapsed, Smeilleur, LeProb->getNomFichier());	//Logguer le temps et la meilleure solution
		dTheBestFitness = 100000;

	}


	return 0;


}
