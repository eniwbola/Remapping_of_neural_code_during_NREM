#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<time.h>
#include<random>
#include "1DConn_fullRandom.h"

using namespace std;

#define step 0.05
// make a scan for gks, change the gks when running the dynamics , so same gks before and after learning but scan over 
// fourier spectrum peak for different hetergeneity size figure 2
// figure 4 do different numbers of inbibitory on, so between a and b, then make figure c for the different grades between
struct neuron
{
	int ID;

	double V_n;
	double m_K;
	double h_Na;
	double s_gKs;

	double gNa, gK, gKs, gLeak;
	double Vrev_Na, Vrev_K, Vrev_leak;

	double I_DC, I_syn_exc, I_syn_inhb;

	bool isInhibitory, isInHet, print, dcOn;

	int printNum, dcCount;

	vector<int> spikeTimes;

	vector<int> connList, preConnList;
	vector<double> connWeights;
	vector<double*> preConnWeights;
};

class printClass
{
private:

	string seedString, gConnString, dcInpString, hetString;
	vector<double> connVals, dcVals, hetVals;

public:
	string getString(int var)
	{
		string returnString;

		switch (var)
		{
		case 0:
			returnString = seedString;
			break;
		case 1:
			returnString = gConnString + dcInpString + hetString;
			break;
		}
		return returnString;
	}

	printClass(int _seed, double _connVals[], double _dcVals[], double _hetVals[], double dcWidth)
	{
		connVals.clear();
		dcVals.clear();

		ostringstream connConv[4], dcConv1, dcConv2, dcConv3, seedConv, hetConv1, hetConv2, dcwConv;
		string connString[4], dcString[3], _hetString[2], dcwString;

		seedConv << _seed;
		seedString = seedConv.str();

		for (int conn_ind = 0; conn_ind < 4; conn_ind++)
		{
			connVals.push_back(_connVals[conn_ind]);
			connConv[conn_ind] << connVals.back();
			connString[conn_ind] = connConv[conn_ind].str();

			unsigned int repConn = connString[conn_ind].find(".");
			if (repConn < connString[conn_ind].size())
			{
				connString[conn_ind].replace(repConn, 1, "p");
			}
		}

		dcConv1 << _dcVals[0];
		dcString[0] = dcConv1.str();
		unsigned int repDC = dcString[0].find(".");
		if (repDC < dcString[0].size())
		{
			dcString[0].replace(repDC, 1, "p");
		}

		dcConv2 << _dcVals[1];
		dcString[1] = dcConv2.str();
		repDC = dcString[1].find(".");
		if (repDC < dcString[1].size())
		{
			dcString[1].replace(repDC, 1, "p");
		}

		dcConv3 << _dcVals[2];
		dcString[2] = dcConv3.str();
		repDC = dcString[2].find(".");
		if (repDC < dcString[2].size())
		{
			dcString[2].replace(repDC, 1, "p");
		}

		double minDC = _dcVals[2] * (1 - dcWidth);
		dcwConv << minDC;

		dcwString = dcwConv.str();
		unsigned int repdcw = dcwString.find(".");
		if (repdcw < dcwString.size())
		{
			dcwString.replace(repdcw, 1, "p");
		}

		gConnString = "EE=" + connString[0] + "_EI=" + connString[1] + "_II=" + connString[2] + "_IE=" + connString[3];
		dcInpString = "_EDC=" + dcString[0] + "_IDC=" + dcString[1] + "_pDC=" + dcString[2] + "to" + dcwString;

		hetConv1 << _hetVals[0];
		hetConv2 << _hetVals[1];

		_hetString[0] = "_hetSize=" + hetConv1.str();
		_hetString[1] = "_hetDC=" + hetConv2.str();
		hetString = _hetString[0] + _hetString[1];

	}
};

double Equations(neuron n, double solns[], int eqn)
{

	double val = 0.;

	switch (eqn)
	{
	case 0:
	{
		double m_inf = 1. / (1. + exp((-solns[0] - 30.) / 9.5));
		double I_Na = -n.gNa*pow(m_inf, 3.)*solns[1] * (solns[0] - n.Vrev_Na);
		double I_K = -n.gK*pow(solns[2], 4.)*(solns[0] - n.Vrev_K);
		double I_Ks = -n.gKs*solns[3] * (solns[0] - n.Vrev_K);
		double I_leak = -n.gLeak*(solns[0] - n.Vrev_leak);

		val = I_Na + I_K + I_leak + I_Ks + n.I_DC - n.I_syn_exc*(solns[0] - 0.) - n.I_syn_inhb*(solns[0] - (-75.));// +_n.I_AC - _n.I_syn_inhb*(solns[0] - (-75.)) - _n.I_syn_exc*(solns[0] - 0.) + _n.I_noise;

		break;
	}
	case 1:
	{
		double h_inf = 1. / (1. + exp((solns[0] + 53.) / 7.));
		double tau_h = 0.37 + 2.78 / (1. + exp((solns[0] + 40.5) / 6.));

		val = (h_inf - solns[1]) / tau_h;
		break;
	}
	case 2:
	{
		double n_inf = 1. / (1. + exp((-solns[0] - 30.) / 10.));
		double tau_n = 0.37 + 1.85 / (1 + exp((solns[0] + 27.) / 15.));

		val = (n_inf - solns[2]) / tau_n;
		break;
	}
	case 3:
	{
		double s_inf = 1. / (1. + exp((-solns[0] - 39.) / 5.));

		val = (s_inf - solns[3]) / 75.;
		break;
	}
	}

	return val;
}

vector<double> synapticPulse(int refracTime, double tau_syn, double tau_spk)
{
	vector<double> sp;
	sp.resize((int)refracTime / step);

	for (int x = 0; x < (int)refracTime / step; x++)
	{
		sp[x] = exp(-((double)x*step - tau_spk) / tau_syn) - exp(-((double)x*step - tau_spk) / 0.3);
	}

	return sp;
}

int main(int argc, char *argv[])
{
	clock_t timeVal = clock();

	int randomSeed = 1750;

	mt19937 mersenne;
	mersenne.seed(randomSeed);

	const int simulationTime = 10e4, transientTime = 1e4, learningTime = 20e4;

	const int networkSize = 1000, eNetworkSize = 800;
	double dcWidth = 1.;
	double connProbs[3] = { 0.3, 0.5, 0.06 };
	double EE = 1e-5, EI = 0.00046, II = 0.0013, IE = 5e-4;//IE= 5e-4;
	double connectivityVals[4] = { EE, EI, II, IE };
	int hetSize = 200;
	double hetStrength = 10.;
	double hetData[2] = { (double)hetSize, hetStrength };
	double learningRate = 1e-5;
	double noiseVal = 0., noiseVec[5] = { 1e-5, 2e-5,5e-5,7e-5,1e-4 };

	double dcVals[3] = { 1.2, -0.085};

	ostringstream seedConv;
	seedConv << randomSeed;

	//string fileBase = "C:/Users/eniwbola/Desktop/Research/Quinton_Proj/noiseSweep/noiseSweep/Data/";  //"D:\\qmskill\\Data\\HH_based\\fastCode\\Het_wNoise\\test\\";
	//string fileBase = "C:/Users/eniwbola/Desktop/Research/Quinton_Proj/sizeSweep/sizeSweep/code/";
	//string fileBase = "/mnt/c/Users/eniwbola/Desktop/Research/Quinton_Proj/sizeSweep/sizeSweep/";
	string fileBase = "/mnt/c/Users/eniwbola/Desktop/Research/Quinton_Proj/noiseSweep/noiseSweep/";
	string paramFileName = fileBase + "ParamData/" + seedConv.str() + ".dat";
    cout << "param_file_name" << paramFileName << endl;
	ofstream paramFileOutput(paramFileName);
	paramFileOutput.close();
	paramFileOutput.open(paramFileName, ios_base::app);

	for (int ns_ind = 0; ns_ind < 5; ns_ind++)
	{
		cout << "ns_ind " << ns_ind << endl; 
		struct neuron *nrn = new neuron[networkSize];

		printClass print(randomSeed, connectivityVals, dcVals, hetData, dcWidth);

		uniform_int_distribution<int> randomNeuron(0, networkSize - 1);
		uniform_real_distribution<double> probability(0., 1.), dist(-1., 1.);

		vector<double> synPulse_long = synapticPulse(1500, 250., 0.08), synPulse_short = synapticPulse(200, 5., 0.);

		vector<string> typeString(3);
		typeString[0] = "_T1NL_";
		typeString[1] = "_T2L_";
		typeString[2] = "_T1AL_";
		string crntType = typeString[0], learningPhase = typeString[1];

		noiseVal = noiseVec[ns_ind];

		paramFileOutput << ns_ind << "\t" << EE << "\t" << EI << "\t" << II << "\t" << IE << "\t" << hetSize << "\t" << hetStrength << "\t" << noiseVal << endl;

		ostringstream nsConv;
		nsConv << ns_ind;

		string spkFile = fileBase + "spikeData/" + seedConv.str() + typeString[0] + nsConv.str() + ".dat";
		string connFile = fileBase + "Connectivity/" + seedConv.str() + typeString[0] + nsConv.str() + ".dat";
		string weightFile_pre = fileBase + "Weights/" + seedConv.str() + typeString[0] + nsConv.str() + ".dat";
		string weightFile_post = fileBase + "Weights/" + seedConv.str() + typeString[2] + nsConv.str() + ".dat";

		/*
		string fileBase = "/scratch/michalz_flux/qmskill/NMDA_Type2_STDP/";
		string spkFile = fileBase + "spikeData/3Part_" + print.getString(0) + "_" + print.getString(1) + ".dat";
		string connFile = fileBase + "Connectivity/3Part_" + print.getString(0) + "_" + print.getString(1) + ".dat";
		string weightFile_pre = fileBase + "Weights/3Part_" + print.getString(0) + "_" + print.getString(1) + ".dat";
		string weightFile_post = fileBase + "Weights/3Part_" + print.getString(0) + "_" + print.getString(1) + ".dat";
		*/
		ofstream conn_output(connFile), spike_output(spkFile), weight_output(weightFile_pre), cont_output[networkSize];

		conn_output.close();
		spike_output.close();
		weight_output.close();

		conn_output.open(connFile, ios_base::app);
		spike_output.open(spkFile, ios_base::app);
		weight_output.open(weightFile_pre, ios_base::app);

		//===============================================================
		//===============================================================
		//Network Initialization

		//int _rSeed, int _netSize, int _eNetSize, int _hetSize, double _rwp, double _connDensity, double *_pVals
		generateConnectivity conn(randomSeed, networkSize, eNetworkSize, 0, 0.2, 0.06, connProbs);

		for (int n_ind = 0; n_ind < networkSize; n_ind++)
		{
			nrn[n_ind].ID = n_ind;

			nrn[n_ind].V_n = -70.*probability(mersenne);// -64.8314;
			nrn[n_ind].h_Na = probability(mersenne);
			nrn[n_ind].m_K = probability(mersenne);
			nrn[n_ind].s_gKs = probability(mersenne);// 0.5;

			nrn[n_ind].gNa = 24.;
			nrn[n_ind].gK = 3.;
			nrn[n_ind].gLeak = 0.02;

			nrn[n_ind].Vrev_Na = 55.;
			nrn[n_ind].Vrev_K = -90.;
			nrn[n_ind].Vrev_leak = -60.;

			if (crntType.compare(1, 1, "2") == 0)
			{
				nrn[n_ind].gKs = 1.5;//0.1 for Type 1, 0.8 for Type 2	
			}
			else
			{
				nrn[n_ind].gKs = 0.;//0.1 for Type 1, 0.8 for Type 2
			}

			nrn[n_ind].I_DC = dcVals[0] - 0.1*dist(mersenne);// (dcVals[0] + 0.2) - 0.4*probability(mersenne);// *(1 - 0.2*probability(mersenne));
			nrn[n_ind].I_syn_exc = 0.;
			nrn[n_ind].I_syn_inhb = 0.;

			nrn[n_ind].preConnList.clear();
			nrn[n_ind].preConnWeights.clear();
			nrn[n_ind].spikeTimes.clear();

			nrn[n_ind].printNum = -1;
			nrn[n_ind].dcCount = -1;

			nrn[n_ind].isInhibitory = conn.isInhibitory[n_ind];
			nrn[n_ind].isInHet = false;

			if (nrn[n_ind].isInhibitory)
			{
				nrn[n_ind].I_DC = dcVals[1] * (1 - 0.1*probability(mersenne));// (dcVals[1] + 0.1) - 0.2*probability(mersenne);
				nrn[n_ind].gKs = 0.;
			}

		}

		int hetCount = 0;

		while (true)
		{
			int randNrn = randomNeuron(mersenne);

			bool choose = false;

			if (!nrn[randNrn].isInhibitory && !nrn[randNrn].isInHet)
			{
				choose = true;
				nrn[randNrn].isInHet = true;
				hetCount++;
			}
			else
			{
				continue;
			}

			if (hetCount == hetSize)
			{
				break;
			}
		}

		for (int n_ind = 0; n_ind < networkSize; n_ind++)
		{

			nrn[n_ind].connList.resize(conn.connections[n_ind].size());
			nrn[n_ind].connWeights.resize(conn.connections[n_ind].size());

			for (int conn_ind = 0; conn_ind < conn.connections[n_ind].size(); conn_ind++)
			{
				int connNode = conn.connections[n_ind][conn_ind];

				nrn[n_ind].connList[conn_ind] = connNode;
				nrn[connNode].preConnList.push_back(n_ind);

				if (nrn[n_ind].isInhibitory & nrn[connNode].isInhibitory)
				{//II
					nrn[n_ind].connWeights[conn_ind] = II;
				}
				else if (nrn[n_ind].isInhibitory)
				{//IE
					nrn[n_ind].connWeights[conn_ind] = IE;
				}
				else if (nrn[connNode].isInhibitory)
				{//EI
					nrn[n_ind].connWeights[conn_ind] = EI;
				}
				else
				{//EE
					nrn[n_ind].connWeights[conn_ind] = EE;

					if (nrn[n_ind].isInHet)
					{
						nrn[n_ind].connWeights[conn_ind] = EE*hetStrength;
					}
				}

				nrn[connNode].preConnWeights.push_back(&nrn[n_ind].connWeights[conn_ind]);
			}
		}

		//Print connections and weights
		for (int n_ind = 0; n_ind < networkSize; n_ind++)
		{

			int nID = n_ind + 1;
			if (nrn[n_ind].isInhibitory)
			{
				nID *= -1;
			}

			conn_output << nID;
			weight_output << nID;

			for (int conn_ind = 0; conn_ind < nrn[n_ind].connList.size(); conn_ind++)
			{
				conn_output << "\t" << nrn[n_ind].connList[conn_ind] + 1;
				weight_output << "\t" << nrn[n_ind].connWeights[conn_ind];
			}

			conn_output << endl;
			weight_output << endl;
		}

		conn_output.close();
		weight_output.close();


		//===============================================================
		//===============================================================
		//Simulation

		mersenne.seed(clock());

		//Loops through time
		for (int t_ind = 0; t_ind <= transientTime + learningTime + 2 * simulationTime; t_ind++)
		{
			vector<int> updateNeurons;
			int modVal = t_ind % (int)20e4;

			//if (modVal == 0 && t_ind > 0 && t_ind < 60e4)
			//If there are phases defined by differing levels of Ach, these statements create new files for data in those phases.
			if (t_ind == simulationTime + transientTime)
			{
				crntType = typeString[1];
				spike_output.close();
				string spkFileName_ = spkFile = fileBase + "spikeData/" + seedConv.str() + typeString[1] + nsConv.str() + ".dat";
				spike_output.open(spkFileName_);
				spike_output.close();
				spike_output.open(spkFileName_, ios_base::app);

				for (int n_ind = 0; n_ind < networkSize; n_ind++)
				{
					if (nrn[n_ind].isInhibitory)
					{
						continue;
					}

					nrn[n_ind].gKs = 1.5;
				}
			}
			else if (t_ind == simulationTime + transientTime + learningTime)
			{
				crntType = typeString[2];
				spike_output.close();
				string spkFileName_ = spkFile = fileBase + "spikeData/" + seedConv.str() + typeString[2] + nsConv.str() + ".dat";
				spike_output.open(spkFileName_);
				spike_output.close();
				spike_output.open(spkFileName_, ios_base::app);

				for (int n_ind = 0; n_ind < networkSize; n_ind++)
				{
					if (nrn[n_ind].isInhibitory)
					{
						continue;
					}

					nrn[n_ind].gKs = 0.;
				}
			}

			//Loops through neurons
			for (int n_ind = 0; n_ind < networkSize; n_ind++)
			{
				nrn[n_ind].dcCount++;

				//Calculate Synaptic Input
				for (int conn_ind = 0; conn_ind < nrn[n_ind].preConnList.size(); conn_ind++)
				{
					int preConnNode = nrn[n_ind].preConnList[conn_ind];

					if (nrn[preConnNode].spikeTimes.size() > 0)
					{
						int timeDiff = t_ind - nrn[preConnNode].spikeTimes.back();

						if (nrn[preConnNode].isInhibitory)
						{

							if (timeDiff <= synPulse_short.size() && timeDiff > 0)
							{
								nrn[n_ind].I_syn_inhb = nrn[n_ind].I_syn_inhb + (*nrn[n_ind].preConnWeights[conn_ind]) * synPulse_short[timeDiff - 1];
							}

						}
						else
						{
							if (nrn[n_ind].isInhibitory)
							{
								if (timeDiff <= synPulse_short.size() && timeDiff > 0)
								{
									nrn[n_ind].I_syn_exc = nrn[n_ind].I_syn_exc + (*nrn[n_ind].preConnWeights[conn_ind]) * synPulse_short[timeDiff - 1];
								}
							}
							else
							{
								if (timeDiff <= synPulse_long.size() && timeDiff > 0)
								{
									nrn[n_ind].I_syn_exc = nrn[n_ind].I_syn_exc + (*nrn[n_ind].preConnWeights[conn_ind]) * synPulse_long[timeDiff - 1];
								}

							}
						}
					}
				}


				//Adds noise to synaptic input. Ideally, this can be added anywhere.
				if (probability(mersenne) <= noiseVal)
				{
					nrn[n_ind].I_syn_exc += 10.;
				}


				//Rk4

				double solutions[4];

				solutions[0] = nrn[n_ind].V_n;
				solutions[1] = nrn[n_ind].h_Na;
				solutions[2] = nrn[n_ind].m_K;
				solutions[3] = nrn[n_ind].s_gKs;

				double k1_0 = Equations(nrn[n_ind], solutions, 0);
				double k1_1 = Equations(nrn[n_ind], solutions, 1);
				double k1_2 = Equations(nrn[n_ind], solutions, 2);
				double k1_3 = Equations(nrn[n_ind], solutions, 3);

				solutions[0] = nrn[n_ind].V_n + 0.5*step*k1_0;
				solutions[1] = nrn[n_ind].h_Na + 0.5*step*k1_1;
				solutions[2] = nrn[n_ind].m_K + 0.5*step*k1_2;
				solutions[3] = nrn[n_ind].s_gKs + 0.5*step*k1_3;

				double k2_0 = Equations(nrn[n_ind], solutions, 0);
				double k2_1 = Equations(nrn[n_ind], solutions, 1);
				double k2_2 = Equations(nrn[n_ind], solutions, 2);
				double k2_3 = Equations(nrn[n_ind], solutions, 3);

				solutions[0] = nrn[n_ind].V_n + 0.5*step*k2_0;
				solutions[1] = nrn[n_ind].h_Na + 0.5*step*k2_1;
				solutions[2] = nrn[n_ind].m_K + 0.5*step*k2_2;
				solutions[3] = nrn[n_ind].s_gKs + 0.5*step*k2_3;

				double k3_0 = Equations(nrn[n_ind], solutions, 0);
				double k3_1 = Equations(nrn[n_ind], solutions, 1);
				double k3_2 = Equations(nrn[n_ind], solutions, 2);
				double k3_3 = Equations(nrn[n_ind], solutions, 3);

				solutions[0] = nrn[n_ind].V_n + step*k3_0;
				solutions[1] = nrn[n_ind].h_Na + step*k3_1;
				solutions[2] = nrn[n_ind].m_K + step*k3_2;
				solutions[3] = nrn[n_ind].s_gKs + step*k3_3;

				double k4_0 = Equations(nrn[n_ind], solutions, 0);
				double k4_1 = Equations(nrn[n_ind], solutions, 1);
				double k4_2 = Equations(nrn[n_ind], solutions, 2);
				double k4_3 = Equations(nrn[n_ind], solutions, 3);

				solutions[0] = nrn[n_ind].V_n + step*(k1_0 + 2.*k2_0 + 2.*k3_0 + k4_0) / 6.;
				solutions[1] = nrn[n_ind].h_Na + step*(k1_1 + 2.*k2_1 + 2.*k3_1 + k4_1) / 6.;
				solutions[2] = nrn[n_ind].m_K + step*(k1_2 + 2.*k2_2 + 2.*k3_2 + k4_2) / 6.;
				solutions[3] = nrn[n_ind].s_gKs + step*(k1_3 + 2.*k2_3 + 2.*k3_3 + k4_3) / 6.;

				//Spike threshold is set at -20 mV
				if (nrn[n_ind].V_n <= -20. && solutions[0] > -20.)
				{
					nrn[n_ind].spikeTimes.push_back(t_ind);

					int nID = n_ind + 1, hetID = 0;

					if (nrn[n_ind].isInhibitory)
					{
						nID *= -1;
					}
					else
					{
						updateNeurons.push_back(n_ind);
					}

					if (nrn[n_ind].isInHet)
					{
						hetID = 1;
					}

					spike_output << nID << "\t" << hetID << "\t" << t_ind*step << endl;// "\t" << nrn[n_ind].I_DC << "\t" << nrn[n_ind].I_syn_exc << "\t" << nrn[n_ind].I_syn_inhb << endl;
				}

				nrn[n_ind].V_n = solutions[0];
				nrn[n_ind].h_Na = solutions[1];
				nrn[n_ind].m_K = solutions[2];
				nrn[n_ind].s_gKs = solutions[3];

				nrn[n_ind].I_syn_exc = 0.;
				nrn[n_ind].I_syn_inhb = 0.;
			}

			//STDP
			if (crntType.compare(learningPhase) == 0)
			{

				for (int n_ind = 0; n_ind < updateNeurons.size(); n_ind++)
				{
					int updateNrn = updateNeurons[n_ind];

					double spikeTime_i = t_ind*step;
					double spikeTime_j = 0.;
					double diffTime_ij = 0.;

					//Loops through the connections pre-synaptic to n_ind and increases the weight
					for (int conn_ind = 0; conn_ind < nrn[updateNrn].preConnList.size(); conn_ind++)
					{
						int preConnNode = nrn[updateNrn].preConnList[conn_ind];

						if (nrn[preConnNode].spikeTimes.size() > 0 && !nrn[preConnNode].isInhibitory)
						{
							spikeTime_j = nrn[preConnNode].spikeTimes.back()*step;
							diffTime_ij = spikeTime_i - spikeTime_j;

							if (diffTime_ij > 1.5 / step)
							{
								*nrn[updateNrn].preConnWeights[conn_ind] = *nrn[updateNrn].preConnWeights[conn_ind] + learningRate*exp(-1.*diffTime_ij / 10.);

								/*if (*nrn[updateNrn].preConnWeights[conn_ind] > 10.*EE && false)
								{
									*nrn[updateNrn].preConnWeights[conn_ind] = 10.*EE;
								}*/
							}
						}
					}

					//Loops through the connections projected by n_ind (what it sends info to) and decreases the weight
					for (int conn_ind = 0; conn_ind < nrn[updateNrn].connList.size(); conn_ind++)
					{
						int connNode = nrn[updateNrn].connList[conn_ind];

						if (nrn[connNode].spikeTimes.size() > 0 && !nrn[connNode].isInhibitory)
						{
							spikeTime_j = nrn[connNode].spikeTimes.back()*step;
							diffTime_ij = spikeTime_i - spikeTime_j;

							if (diffTime_ij > 1.5 / step)
							{
								nrn[updateNrn].connWeights[conn_ind] = nrn[updateNrn].connWeights[conn_ind] - learningRate*exp(-1.*diffTime_ij / 10.);

								if (nrn[updateNrn].connWeights[conn_ind] < 0.)
								{
									nrn[updateNrn].connWeights[conn_ind] = 0.;
								}
							}
						}
					}
				}

			}

			updateNeurons.clear();

		}

		//prints final weights to the After Learning (AL) file

		weight_output.open(weightFile_post);

		for (int n_ind = 0; n_ind < networkSize; n_ind++)
		{
			weight_output << n_ind;

			for (int conn_ind = 0; conn_ind < nrn[n_ind].connList.size(); conn_ind++)
			{
				weight_output << "\t" << nrn[n_ind].connWeights[conn_ind];
			}

			weight_output << endl;
		}

		weight_output.close();
		spike_output.close();

		delete[]nrn;

		cout << "Job done..." << endl;
		cout << "Job took " << (double)(clock() - timeVal) / (double)CLOCKS_PER_SEC << " seconds..." << endl;

	}
}
