#include<vector>
#include<iostream>
#include<random>
#include<time.h>

using namespace std;

class generateConnectivity
{
private:
	int rSeed, netSize, eNetSize, iNetSize, hetSize;
	double rwp, connDensity, *connProbs;

	void createConnectivity();
	int findMin(int nInd, vector<int> pop);

public:

	vector<vector<int> > connections, EIPopulations;
	vector<bool> isInhibitory, isInHet;

	int PopID(int testID);

	generateConnectivity(int _rSeed, int _netSize, int _eNetSize, int _hetSize, double _rwp, double _connDensity, double *_pVals)
		:rSeed(_rSeed), netSize(_netSize), eNetSize(_eNetSize), iNetSize(_netSize - _eNetSize), hetSize(_hetSize), rwp(_rwp), connDensity(_connDensity), connProbs(_pVals)
	{
		connections.clear();
		EIPopulations.clear();
		isInhibitory.clear();
		isInHet.clear();

		createConnectivity();
	}
};

int generateConnectivity::findMin(int nInd, vector<int> pop)
{
	int min_ID_diff = netSize*netSize, minID = 0;

	for (int ind = 0; ind < pop.size(); ind++)
	{
		int nodeID = pop[ind];
		if (abs(nInd - nodeID) < min_ID_diff)
		{
			min_ID_diff = abs(nInd - nodeID);
			minID = ind;
		}
	}

	return minID;
}

int generateConnectivity::PopID(int testID)
{
	int _popID_ = -1;

	for (int n_ind = 0; n_ind < EIPopulations[1].size(); n_ind++)
	{
		if (testID == EIPopulations[1][n_ind])
		{
			_popID_ = 1;
			return _popID_;
		}
	}

	for (int n_ind = 0; n_ind < EIPopulations[0].size(); n_ind++)
	{
		if (testID == EIPopulations[0][n_ind])
		{
			_popID_ = 0;
			return _popID_;
		}
	}

	return _popID_;
}

void generateConnectivity::createConnectivity()
{
	mt19937 mers(rSeed + 1984);

	uniform_real_distribution<double> probability(0., 1.);
	uniform_int_distribution<int> randomNode(0, netSize - 1);

	vector< vector<int> > rewire;
	vector<int> inhibitoryPopulation, excitatoryPopulation;

	int inhbVar = netSize / iNetSize;

	isInhibitory.resize(netSize);
	isInHet.resize(netSize);

	//Identifies inhibitory neurons and neurons belonging to heterogeneity
	for (int n_ind = 0; n_ind < netSize; n_ind++)
	{
		if (n_ind % inhbVar == 0)
		{//Places inhibitory neurons evenly throughout the network
			isInhibitory[n_ind] = true;
			inhibitoryPopulation.push_back(n_ind);
		}
		else
		{
			isInhibitory[n_ind] = false;
			excitatoryPopulation.push_back(n_ind);
		}
	}

	EIPopulations.push_back(excitatoryPopulation);
	EIPopulations.push_back(inhibitoryPopulation);

	uniform_int_distribution<int> rand_exc_node(0, excitatoryPopulation.size() - 1);
	uniform_int_distribution<int> rand_inhib_node(0, inhibitoryPopulation.size() - 1);

	//Initial Connectivity Setup
	for (int n_ind = 0; n_ind < netSize; n_ind++)
	{
		//Temporary vectors
		vector<int> conns;
		vector<int> rewire_i;

		if (isInhibitory[n_ind])
		{//Inhibitory connections

			int connNode = 0;
			double val = 0.;

			for (int conn_ind = 0; conn_ind < excitatoryPopulation.size(); conn_ind++)
			{//Loops through the network --- only for I-X connections
				
				connNode = excitatoryPopulation[conn_ind];

				if (probability(mers) <= connProbs[0])
				{
					conns.push_back(connNode);
				}

			}

			for (int conn_ind = 0; conn_ind < inhibitoryPopulation.size(); conn_ind++)
			{
				connNode = inhibitoryPopulation[conn_ind];

				if (connNode == n_ind)
				{
					continue;
				}

				if (probability(mers) <= connProbs[1])
				{
					conns.push_back(connNode);
				}
			}

		}
		else
		{//excitatory connections

			for (int conn_ind = 0; conn_ind < excitatoryPopulation.size(); conn_ind++)
			{
				int connNode = excitatoryPopulation[conn_ind];

				if (n_ind == connNode)
				{
					continue;
				}
				else if (probability(mers) <= connProbs[2])
				{
					conns.push_back(connNode);

					if (probability(mers) <= rwp)
					{
						rewire_i.push_back(conns.size() - 1);
					}
				}
			}

			for (int conn_ind = 0; conn_ind < inhibitoryPopulation.size(); conn_ind++)
			{
				int connNode = inhibitoryPopulation[conn_ind];

				if (probability(mers) <= connProbs[2])
				{
					conns.push_back(connNode);

					if (probability(mers) <= rwp)
					{
						rewire_i.push_back(conns.size() - 1);
					}
				}
			}
		}

		//Adds the temporary connections and rewiring to the master lists
		connections.push_back(conns);
		rewire.push_back(rewire_i);

		conns.clear();
		rewire_i.clear();
	}

	
	//Rewires connections flagged during initial connectivity setup
	//for (int n_ind = 0; n_ind < netSize; n_ind++)
	//{
	//	if (rewire[n_ind].size() == 0)
	//	{//Skips connectivity lists that have no connections to be rewired
	//		continue;
	//	}

	//	for (int ind = 0; ind < rewire[n_ind].size(); ind++)
	//	{//Loops through the connections to be rewired

	//		//Identifies the original node ID
	//		int ogNode = connections[n_ind][rewire[n_ind][ind]];

	//		while (true)
	//		{//Break on a succesful reconnection
	//			bool choose = true;
	//			int chosenNode = 0;

	//			if (ogNode % inhbVar == 0)
	//			{//Rewire to another inhibitory neuron

	//				chosenNode = inhibitoryPopulation[rand_inhib_node(mers)];

	//				for (int conn_ind = 0; conn_ind < connections[n_ind].size(); conn_ind++)
	//				{//Loops through the connections to find duplicate connections
	//					if (chosenNode == connections[n_ind][conn_ind])
	//					{//rejects the proposed new connection
	//						choose = false;
	//						break;
	//					}
	//				}
	//			}
	//			else
	//			{//Rewire to another excitatory neuron
	//			
	//				chosenNode = excitatoryPopulation[rand_exc_node(mers)];

	//				for (int conn_ind = 0; conn_ind < connections[n_ind].size(); conn_ind++)
	//				{//Lopps through the connections to find duplicate connections
	//					if (chosenNode == connections[n_ind][conn_ind])
	//					{//rejects the proposed new connection
	//						choose = false;
	//						break;
	//					}
	//				}
	//			}

	//			if (choose)
	//			{//Rewires the connection and break the continuous loop
	//				connections[n_ind][ogNode] = chosenNode;
	//				break;
	//			}
	//		}
	//	}
	//}

}
