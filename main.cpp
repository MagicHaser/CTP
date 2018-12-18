#include <iostream>
#include <ilcplex/ilocplex.h>
#include <iomanip>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include<vector>
#include <string>



#include <fstream>
using namespace std;

#define RC_EPS 1.0e-6

int max_pool_size = 196;
double actualUtility_best = 0.0; 
int cut_count = 0;
int userCut_count = 0;
double*** bestSolVol = NULL;
double*** bestSolVolOtherSetting = NULL;


struct milktype {
	int cal_lvl;
	int bottle_size;
	//double target;
	double cal_lb;
	double cal_ub;
	//double coe_surplus;
	//int target_range;
	double protein_lb;
	int productionTarget;
	double unitUtility_def;
	double unitUtility_sur;
};

struct deposit{
	double volume;
	bool new_mom;
	double calorie;
	double protein;
	double bacteriaFreeProb;
};

struct pool {
	bool is_used;
	int type_idx;
	int volume;
	vector<int> depo_idx;
	double clean_prob;
};



typedef IloArray<IloNumVarArray> IloNumVar2dArray;
typedef IloArray<IloNumVar2dArray> IloNumVar3dArray;
typedef IloArray<IloArray<IloExpr>> IloExpr2dArray;

//IloBool OptCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const IloNumArray depoProb, const IloNumArray typeUtility, IloCplex cplex, IloExpr v_singlecut, IloExpr2dArray v_multicut, IloNum& singlecutRHS, IloNumArray multicutRHS);
IloBool OptSingleCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const vector<deposit> depo, const vector<milktype> type, IloCplex cplex, IloExpr v_singlecut, IloNum& singlecutRHS, vector<pool> pools);
IloBool OptMultiCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const vector<deposit> depo, const vector<milktype> type, IloCplex cplex, IloArray<IloExpr> v_multicut, IloNumArray multicutRHS, vector<pool> pools);
IloBool OptAdditionalCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const vector<deposit> depo, const vector<milktype> type, IloCplex cplex, IloArray<IloExpr> v_multicut, IloNum& multicutRHS, vector<pool> pools, IloBool & bestUtilityImproved);
IloNumArray meetTargetProb(const IloNumArray3 allocateVolVarsSol, vector<deposit> depo, vector<milktype> type);
void ObtainPoolInfo(const IloNumArray3 allocateVolVarsSol, vector<deposit> depo, vector<pool> &pools);
//void UpdatePoolProb(vector<pool> &pools, IloNumArray depo_prob);
//IloBoolArray ObtainTypeProduction(IloNumArray2 &TypePoolVol, IloNumArray2 &TypePoolProb, vector<pool> pools);
IloNumArray CalculateMeetTargetProb(vector<pool> pools, IloNumArray productionTargets, IloNumArray depoCleanProb, IloNumArray depoVol);
IloNumArray WhetherMeetTarget(vector<pool> pools, IloNumArray productionTargets, IloBoolArray depo_clean) ;
IloNumArray SumofTwoArray(IloNumArray array1, IloNum prob1, IloNumArray array2, IloNum prob2);
IloNum evaluateActualUtility(IloNumArray depoCleanProb, vector<pool> pools, vector<milktype> type);
double evaluateUtilityBySol(vector<milktype> type, vector<double> depoCleanProb, double *** Sol);
//IloNum ObtainMeetTargetProb(IloNumArray poolVol, IloNumArray poolProb, IloNum productionTarget);
//IloNumArray ObtainMeetTargetProbAll(IloNumArray2 TypePoolVol, IloNumArray2 TypePoolProb, IloNumArray productionTarget, IloBoolArray producedType);


ILOLAZYCONSTRAINTCALLBACK7(LShapedLazyCallback, IloNumVar3dArray, allocateVolVars, IloNumVar2dArray, allocateBinVars, IloNumVarArray, expectedUtilityType, vector<deposit>, depo, vector<milktype>, type, vector<pool>, pools, IloCplex, cplex)
{
	IloInt i, j, n;
	IloEnv masterEnv = getEnv();
	IloInt num_pool = allocateVolVars[0].getSize();
	IloInt num_depo = allocateVolVars.getSize();
	IloInt num_type = expectedUtilityType.getSize() - 1;
	IloInt cutType = 4; //1: single cut; 2: multi cuts by type; 3: duplicate cuts by enumerate 



	// Get the current x solution

	IloNumArray3 allocateVolVarsSol(masterEnv, num_depo);
	for (i = 0; i < num_depo; i++) {
		allocateVolVarsSol[i] = IloNumArray2(masterEnv, num_pool);
		for (j = 0; j < num_pool; j++) {
			allocateVolVarsSol[i][j] = IloNumArray(masterEnv);
			getValues(allocateVolVarsSol[i][j], allocateVolVars[i][j]);
		}
	}

	IloNumArray2 allocateBinVarsSol(masterEnv, num_depo);
	for (i = 0; i < num_depo; i++) {
		allocateBinVarsSol[i] = IloNumArray(masterEnv);
		getValues(allocateBinVarsSol[i], allocateBinVars[i]);
	}

	IloNumArray expectedUtilityTypeSol(masterEnv);
	getValues(expectedUtilityTypeSol, expectedUtilityType);

	// Optimality cut


//	IloBool cutStat = OptCutAdded(allocateVolVars, allocateVolVarsSol, allocateBinVars, allocateBinVarsSol, expectedUtilityType, expectedUtilityTypeSol, max_pool_size, depoProb, typeUtility, cplex, v_singlecut, v_multicut, singlecutRHS, multicutRHS);

	IloBool cutStat;
	if (cutType == 1 || cutType == 5) {
		IloExpr v_singlecut(masterEnv);
		IloNum singlecutRHS;
		cutStat = OptSingleCutAdded(allocateVolVars, allocateVolVarsSol, allocateBinVars, allocateBinVarsSol, expectedUtilityType, expectedUtilityTypeSol, max_pool_size, depo, type, cplex, v_singlecut, singlecutRHS, pools);
		if (cutStat) {
			//masterEnv.out()<< v_singlecut << " >= " <<  singlecutRHS <<endl; 
			add(v_singlecut >= singlecutRHS).end();
			//add(expectedUtilityType[num_type] >= actualUtility_best).end();
			/*		for (j = 0; j < num_pool; j++) {
			for (j1 = 0; j1<num_pool; j1++)
			add(v_multicut[j][j1]>= multicutRHS[j]).end();
			}	*/
		}

		v_singlecut.end();
		if (cutType == 5) {
			IloArray<IloExpr> v_multicut(masterEnv, num_type);
			IloNumArray multicutRHS(masterEnv, num_type);
			for (n = 0; n < num_type; n++) {
				v_multicut[n] = IloExpr(masterEnv);
			}
			cutStat = OptMultiCutAdded(allocateVolVars, allocateVolVarsSol, allocateBinVars, allocateBinVarsSol, expectedUtilityType, expectedUtilityTypeSol, max_pool_size, depo, type, cplex, v_multicut, multicutRHS, pools);
			if (cutStat) {
				for (n = 0; n < num_type; n++)
					add(v_multicut[n] >= multicutRHS[n]).end();
				//add(expectedUtilityType[num_type] >= actualUtility_best).end();
			}
			for (n = 0; n < num_type; n++) {
				v_multicut[n].end();
			}
			v_multicut.end();
			multicutRHS.end();
		}
	}
	else if (cutType == 2) {
		IloArray<IloExpr> v_multicut(masterEnv, num_type);
		IloNumArray multicutRHS(masterEnv, num_type);
		for (n = 0; n < num_type; n++) {
			v_multicut[n] = IloExpr(masterEnv);
		}
		cutStat = OptMultiCutAdded(allocateVolVars, allocateVolVarsSol, allocateBinVars, allocateBinVarsSol, expectedUtilityType, expectedUtilityTypeSol, max_pool_size, depo, type, cplex, v_multicut, multicutRHS, pools);
		if (cutStat) {
			for (n = 0; n < num_type; n++)
				add(v_multicut[n] >= multicutRHS[n]).end();
			add(expectedUtilityType[num_type] >= actualUtility_best).end();
		}
		for (n = 0; n < num_type; n++) {
			v_multicut[n].end();
		}
		v_multicut.end();
		multicutRHS.end();
	}
	else if (cutType == 3 || cutType == 4) {
		IloArray<IloExpr> v_additioncut(masterEnv, num_pool);
		IloNum additioncutRHS;
		for (j = 0; j < num_pool; j++) {
			v_additioncut[j] = IloExpr(masterEnv);
		}
		IloBool bestUtilityImproved = IloFalse;
		cutStat = OptAdditionalCutAdded(allocateVolVars, allocateVolVarsSol, allocateBinVars, allocateBinVarsSol, expectedUtilityType, expectedUtilityTypeSol, max_pool_size, depo, type, cplex, v_additioncut, additioncutRHS, pools, bestUtilityImproved);
		if (bestUtilityImproved)
			add(expectedUtilityType[num_type] >= actualUtility_best).end();
		if (cutStat) {
			for (j = 0; j < num_pool; j++)
				add(v_additioncut[j] >= additioncutRHS).end();
		}
		for (j = 0; j < num_pool; j++) {
			v_additioncut[j].end();
		}
		v_additioncut.end();
		if (cutType == 4) {
			IloArray<IloExpr> v_multicut(masterEnv, num_type);
			IloNumArray multicutRHS(masterEnv, num_type);
			for (n = 0; n < num_type; n++) {
				v_multicut[n] = IloExpr(masterEnv);
			}
			cutStat = OptMultiCutAdded(allocateVolVars, allocateVolVarsSol, allocateBinVars, allocateBinVarsSol, expectedUtilityType, expectedUtilityTypeSol, max_pool_size, depo, type, cplex, v_multicut, multicutRHS, pools);
			if (cutStat) {
				for (n = 0; n < num_type; n++)
					add(v_multicut[n] >= multicutRHS[n]).end();
				//add(expectedUtilityType[num_type] >= actualUtility_best).end();
			}
			for (n = 0; n < num_type; n++) {
				v_multicut[n].end();
			}
			v_multicut.end();
			multicutRHS.end();
		}
	}
	else {
		masterEnv.out() << "Unspecific cut code" << endl;
	}

	for (i = 0; i < num_depo; i++) {
		for (j = 0; j < num_pool; j++) {
			allocateVolVarsSol[i][j].end();
		}
		allocateVolVarsSol[i].end();
	}
	allocateVolVarsSol.end();

	for (i = 0; i < num_pool; ++i) {
		allocateBinVarsSol[i].end();
	}
	allocateBinVarsSol.end();
	expectedUtilityTypeSol.end();

	return;

} // END LShapedLazyCallback


int main(void)
{
	string filefolder;
	filefolder = "C:\\user\\rus19\\Dropbox\\Dropbox\\Pitt\\Research\\MilkBank\\dataset\\Prepooling\\";

	int scale = 100;
	double min_portion = 600 / scale;

	int max_split, max_depoperpool;
	double unitUtility_20 = 1;
	double lexi_ratio = 3;//3
	double AddUtility = 2;//2

	double penaltySplit = 0.01;
	double cost_pool = 1;
	bool additional_cut = 1;
	bool anti_symm_vol = 0;
	bool anti_symm_depoperpool = 0;
	bool anti_symm_hierarchy = 1;
	bool minimize_num_pool = 0;
	bool multi_cut = 1;
	bool set_priority_u = 1;
	bool set_priority_x = 0;
	bool evaluate_keptalone = 1;
	bool minimize_split = 1;
	bool separateFracSols = 0;

	ifstream inputfile;
	ofstream LogFile(filefolder + "log8_new.txt");
	ofstream outputfile(filefolder + "result8_new.txt");
	ofstream resultfile(filefolder + "result8_new.csv");

	for (int input_idx = 1; input_idx < 4; input_idx++){

	vector<deposit> depo;
	vector<milktype> type;
	vector<pool> pools;


				int num_pool;

	double volume, protein, calorie, probability;
	int type_cal, type_bot_size, production_target;

	bool new_mom;
	double type_cal_lb, type_cal_ub, type_protein_lb, unitUtility_def, unitUtility_sur;
	char type_beg, type_end, next_char;

		char input_idx_char = '0' + input_idx;
		int num_established = 4;
		int num_newmom = 3;
		for (int mm = 1; mm < num_established+1; mm++){//5
			for (int m = 1; m < num_newmom+1; m++){//4
				actualUtility_best = 0.0;
				cut_count = 0;
				userCut_count = 0;

				cout << mm << " mm:m " << m << endl;
				outputfile << mm << " mm:m " << m << endl;
				double established_mom_prob = 0.025*mm;
				double new_mom_prob = 0.05*mm*m;

				resultfile << established_mom_prob << "," << new_mom_prob << ",";

				inputfile.open(filefolder + "input" + input_idx_char + ".txt");
				//char input_idx_char[1];
				//itoa(input_idx, input_idx_char, 10);
				//if (!numerical_study || input_idx == 0) inputfile.open(txtfilefolder + "input.txt");
				//else inputfile.open(txtfilefolder + "input" + input_idx_char + ".txt");	inputfile.open("C:/user/rus19/Dropbox/Dropbox/Pitt/Research/MilkBank/dataset/TestRun/Week0829_offline/Solver_input/08302016_modobj.txt");

				if (!inputfile.is_open()) {
					cout << "Cannot generate optimization problem. Please check the input.\n";
					return -1;
				}
				else {
					string word;
					char buffer_alt[10];

					inputfile >> buffer_alt;
					outputfile << " --------------- " << buffer_alt << " --------------- " << endl;
					inputfile >> type_beg;
					do
					{
						type_protein_lb = 0;
						inputfile >> type_cal >> type_bot_size >> production_target >> type_cal_lb >> type_cal_ub >> type_protein_lb >> type_end >> type_beg;
						unitUtility_def = (1 + (type_cal == 20)*(AddUtility - 1))*unitUtility_20*(lexi_ratio);
						unitUtility_sur = unitUtility_20 / (1 + (type_cal != 20)*(lexi_ratio - 1));
						//cout << type_cal << "  " << type_bot_size << endl;
						milktype next_type = { type_cal, type_bot_size, type_cal_lb, type_cal_ub, type_protein_lb, production_target*type_bot_size/100, unitUtility_def, unitUtility_sur };
						type.push_back(next_type);
					} while (type_beg == '{');
					max_depoperpool = type_beg - '0';
					inputfile >> max_split;

					while (inputfile >> volume >> next_char >> calorie >> protein)
					{
						//cout << next_char << endl;
						//if(protein_rqmt) protein = 0;
						if (next_char == 'y' || next_char == 'Y') {
							probability = new_mom_prob;
							new_mom = true;
						}
						else {
							probability = established_mom_prob;
							new_mom = false;
							//cout << "vol: "<< volume  <<" " << calorie  << endl;
						}
						deposit next = { volume / scale, new_mom, calorie, protein, 1 - probability };
						depo.push_back(next);
						//inputfile >> next_char;
					}
					inputfile.close();
					num_pool = depo.size();
					for (int j = 0; j < num_pool; j++) {
						vector<int> depo_idx_new;
						pool nextPool = { false, -1, 0, depo_idx_new, 1.0 };
						pools.push_back(nextPool);
					}

				}



				IloEnv   env;

				try {
					IloModel model(env), model_keptalone(env);
					IloCplex cplex(env), cplex_keptalone(env);

					IloInt num_depo = depo.size();

					IloInt num_type = type.size();
					IloInt i, j, n, n1;
					//		for (i = 0; i < num_depo; i++){
					//			bestSolVol[i] = IloNumArray2(env, num_pool);
					//			for (j = 0; j < num_pool; j++) 
					//				bestSolVol[i][j] = IloNumArray(env, num_type);
					//		}



					bestSolVol = new double**[(int)num_depo];
					for (i = 0; i < num_depo; i++) {
						bestSolVol[i] = new double*[(int)num_pool];
						for (j = 0; j < num_pool; j++) {
							bestSolVol[i][j] = new double[(int)num_type];
							for (n = 0; n < num_type; n++) {
								bestSolVol[i][j][n] = 0;
							}
						}
					}

					if (m*mm == 1) {
						bestSolVolOtherSetting = new double**[(int)num_depo];
						for (i = 0; i < num_depo; i++) {
							bestSolVolOtherSetting[i] = new double*[(int)num_pool];
							for (j = 0; j < num_pool; j++) {
								bestSolVolOtherSetting[i][j] = new double[(int)num_type];
								for (n = 0; n < num_type; n++) {
									bestSolVolOtherSetting[i][j][n] = 0;
								}
							}
						}
					}


					/*
					IloNumVarArray capVars(env, nbWhouses, 0, 10, ILOINT); // Used capacities
					IloNumArray    capLbs(env, nbWhouses, 2, 3, 5, 7); // Minimum usage level
					IloNumArray    costs(env, nbWhouses, 1, 2, 4, 6); // Cost per warehouse

					// These variables represent the assigninment of a
					*/														  // load to a warehouse.
					IloNumVar3dArray allocateVolVars(env, num_depo);
					IloNumVar2dArray allocateBinVars(env, num_depo);
					IloNumVar2dArray pooltypeBinVars(env, num_pool);
					IloNumVar2dArray productionByType(env, 2);
					//IloNumVar expectedUtility(env, -IloInfinity, IloInfinity, "psi");
					IloNumVarArray expectedUtilityType(env, num_type + 1, 0, IloInfinity);
					IloNumArray depoProb(env, num_depo);
					IloNumArray typeUtility(env, num_type);


					//define variables
					for (i = 0; i < num_depo; i++) {
						allocateVolVars[i] = IloNumVar2dArray(env, num_pool);
						for (j = 0; j < num_pool; j++) {
							allocateVolVars[i][j] = IloNumVarArray(env, num_type, 0, max_pool_size, ILOINT);
							for (n = 0; n < num_type; n++) {
								char varName[100];
								sprintf_s(varName, "w_%d_%d_%d", (int)i + 1, (int)j + 1, (int)n + 1);
								allocateVolVars[i][j][n].setName(varName);
							}
						}
					}

					for (i = 0; i < num_depo; i++) {
						allocateBinVars[i] = IloNumVarArray(env, num_pool, 0, 1, ILOINT);
						for (j = 0; j < num_pool; j++) {
							char varName[100];
							sprintf_s(varName, "u_%d_%d", (int)i + 1, (int)j + 1);
							allocateBinVars[i][j].setName(varName);
						}
					}

					for (j = 0; j < num_pool; j++) {
						pooltypeBinVars[j] = IloNumVarArray(env, num_type, 0, 1, ILOINT);
						for (n = 0; n < num_type; n++) {
							char varName[100];
							sprintf_s(varName, "x_%d_%d", (int)j + 1, (int)n + 1);
							pooltypeBinVars[j][n].setName(varName);
						}
					}
					for (n = 0; n < num_type; n++) {
						char varName[10];
						sprintf_s(varName, "psi_%d", (int)n + 1);
						expectedUtilityType[n].setName(varName);
					}
					expectedUtilityType[num_type].setName("psi");

					productionByType[0] = IloNumVarArray(env, num_type, 0, IloInfinity);
					productionByType[1] = IloNumVarArray(env, num_type, 0, IloInfinity);

					//env.out() << "after defining variables" << endl;

					//env.out() << allocateVolVars[0].getSize() << allocateVolVars.getSize() <<endl;

					if (multi_cut) {
						IloExpr v(env);
						for (n = 0; n < num_type; n++)
							v += expectedUtilityType[n];
						model.add(v == expectedUtilityType[num_type]);
						v.end();
					}

					// Constraint: calorie and protein requirement for each type
					//			vectorSumConst(depo_Cal, -type[n].lb);
					for (j = 0; j < num_pool; j++) {
						for (n = 0; n < num_type; n++) {
							IloExpr v1(env), v2(env), v3(env);
							for (i = 0; i < num_depo; i++){
								v1 += (depo[i].calorie - type[n].cal_lb) * allocateVolVars[i][j][n];
								v2 += (depo[i].calorie - type[n].cal_ub) * allocateVolVars[i][j][n];
								v3 += (depo[i].protein - type[n].protein_lb) * allocateVolVars[i][j][n];
								//env.out() << depo[i].calorie << " " << type[n].cal_lb << " " <<depo[i].calorie - type[n].cal_lb<<endl;
							}
							model.add(v1 >= 0);
							model.add(v2 <= 0);
							model.add(v3 >= 0);
							v1.end();
							v2.end();
							v3.end();
						}
					}
					//env.out() << "reach here" << endl;
					// Constraint: one type per pool
					for (j = 0; j < num_pool; j++) {
						model.add(IloSum(pooltypeBinVars[j]) <= 1);
					}

					// Constraint: bound the pool size
					for (j = 0; j < num_pool; j++) {
						for (n = 0; n < num_type; n++) {
							IloExpr v(env);
							for (i = 0; i < num_depo; i++) {
								v += allocateVolVars[i][j][n];
							}
							model.add(v <= max_pool_size * pooltypeBinVars[j][n]);
							model.add(v >= min_portion * pooltypeBinVars[j][n]);
							v.end();
						}
					}

					// Constraint: limit minimum portion size
					for (i = 0; i < num_depo; i++) {
						for (j = 0; j < num_pool; j++) {
							model.add(IloSum(allocateVolVars[i][j]) >= min_portion * allocateBinVars[i][j]);
							model.add(IloSum(allocateVolVars[i][j]) <= depo[i].volume* allocateBinVars[i][j]);
						}
					}

					// Constraint: process all thawed milk
					for (i = 0; i < num_depo; i++) {
						IloExpr v(env);
						for (j = 0; j < num_pool; j++)
							for (n = 0; n < num_type; n++)
								v += allocateVolVars[i][j][n];
						model.add(v == depo[i].volume);
						v.end();
					}

					// Constraint: limit # of portions to be split into
					for (i = 0; i < num_depo; i++) {
						model.add(IloSum(allocateBinVars[i]) <= max_split);
					}

					// Constraint: limit # of deposit/pool
					for (j = 0; j < num_pool; j++) {
						IloExpr v(env);
						for (i = 0; i < num_depo; i++)
							v += allocateBinVars[i][j];
						model.add(v <= max_depoperpool);
						v.end();
					}

					//Constraint: anti-symmetry
					if (anti_symm_vol){
						for (j = 1; j < num_pool; j++) {
							for (n = 0; n < num_type; n++) {
								IloExpr v(env);
								for (n1 = 0; n1 < n + 1; n1++)
									v += pooltypeBinVars[j - 1][n1];
								model.add(v >= pooltypeBinVars[j][n]);
								v.end();
							}
						}

						for (n = 0; n < num_type; n++) {
							for (j = 1; j < num_pool; j++) {
								IloExpr v(env);
								for (i = 0; i < num_depo; i++) {
									v += allocateVolVars[i][j][n] - allocateVolVars[i][j - 1][n];
								}
								v += max_pool_size * pooltypeBinVars[j - 1][n];
								model.add(v <= max_pool_size);
								v.end();
							}
						}
					}
					if (anti_symm_depoperpool){
						for (j = 1; j < num_pool; j++) {
							IloExpr v(env);
							for (i = 0; i < num_depo; i++) {
								v += allocateBinVars[i][j] - allocateBinVars[i][j - 1];
							}
							model.add(v <= 0);
							v.end();
						}
					}
					if (anti_symm_hierarchy) {
						for (j = 1; j < num_pool; j++) {
							IloExpr v(env);
							for (i = 0; i < num_depo; i++) {
								int ii = i;
								v += pow(2.0, ii + 1) * allocateBinVars[i][j - 1] - pow(2.0, ii + 1) * allocateBinVars[i][j];
							}
							//model.add(v >= 0);
							for (n = 0; n < num_type; n++) {
								model.add(v - pooltypeBinVars[j-1][n] - pooltypeBinVars[j][n] >= -1);
							}


							v.end();
						}

					}


					IloNum actualUtility = 0.0;
					IloNum sumUnitUtility = 0.0;
					int iteration = 0;

					IloNumArray3 allocateVolVarsSol_keptalone(env, num_depo);
					if (evaluate_keptalone) {
						model_keptalone.add(model);
						//model_keptalone.remove(obj);
						for (j = 0; j < num_pool; j++) {
							IloExpr v(env);
							for (i = 0; i < num_depo; i++)
								v += (depo[i].new_mom ? max_depoperpool : 1)* allocateBinVars[i][j];
							model_keptalone.add(v <= max_depoperpool);
							v.end();
						}

						for (n = 0; n < num_type; n++) {
							IloExpr v(env);
							for (i = 0; i < num_depo; i++)
								for (j = 0; j < num_pool; j++)
									v += allocateVolVars[i][j][n];
							v -= productionByType[0][n];
							v -= productionByType[1][n];
							model_keptalone.add(v == 0);
							model_keptalone.add(productionByType[0][n] <= type[n].productionTarget);
							v.end();
						}


						IloExpr obj_kept(env);
						for (n = 0; n < num_type; n++)
							obj_kept += type[n].unitUtility_def*productionByType[0][n] + type[n].unitUtility_sur*productionByType[1][n];
						if (minimize_split)
							for (i = 0; i < num_depo; i++)
								obj_kept -= penaltySplit * IloSum(allocateBinVars[i]);
						model_keptalone.add(IloMaximize(env, obj_kept));
						model.remove(obj_kept);
						obj_kept.end();

						cplex_keptalone.extract(model_keptalone);
						cplex_keptalone.exportModel((filefolder + "model_keptalone.lp").c_str());
						cplex_keptalone.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.001);
						cplex_keptalone.setOut(env.getNullStream());//setout
						cplex_keptalone.setWarning(env.getNullStream());
						cplex_keptalone.setParam(IloCplex::TiLim, 600);
						cplex_keptalone.solve();

						for (n = 0; n < num_type; n++){
							env.out() << "Type " << type[n].cal_lvl << " target is " << type[n].productionTarget << " production is " << cplex_keptalone.getValue(productionByType[0][n]) + cplex_keptalone.getValue(productionByType[1][n]) << " total utility is " << cplex_keptalone.getObjValue() << endl;
							outputfile << "Type " << type[n].cal_lvl << " target is " << type[n].productionTarget << " production is " << cplex_keptalone.getValue(productionByType[0][n]) + cplex_keptalone.getValue(productionByType[1][n]) << " total utility is " << cplex_keptalone.getObjValue() << endl;
						}
						//IloNumArray3 allocateVolVarsSol(env, num_depo);
						double *** solKeptalone = new double **[num_depo];
						for (i = 0; i < num_depo; i++) {
							allocateVolVarsSol_keptalone[i] = IloNumArray2(env, num_pool);
							solKeptalone[i] = new double *[num_pool];
							for (j = 0; j < num_pool; j++) {
								solKeptalone[i][j] = new double [num_type];
								allocateVolVarsSol_keptalone[i][j] = IloNumArray(env);
								cplex_keptalone.getValues(allocateVolVars[i][j], allocateVolVarsSol_keptalone[i][j]);
								for (n = 0; n < num_type; n++) {
									solKeptalone[i][j][n] = cplex_keptalone.getValue(allocateVolVars[i][j][n]);
								}
							}
						}
						//env.out << "after extracting " << endl;

						ObtainPoolInfo(allocateVolVarsSol_keptalone, depo, pools);
						IloNumArray depoCleanProb(env, num_depo);
						vector<double> depoCleanProb_keptalone;
						for (int i = 0; i < num_depo; i++){
							depoCleanProb[i] = depo[i].bacteriaFreeProb;
							depoCleanProb_keptalone.push_back(depo[i].bacteriaFreeProb);
						}
						//env.out() << "before actual" << endl;
						IloNum heu_actualUtility = evaluateActualUtility(depoCleanProb, pools, type);
						//env.out() << "after actual" << endl;
						env.out() << heu_actualUtility <<" kept-alone solution " << evaluateUtilityBySol(type, depoCleanProb_keptalone, solKeptalone)  << endl;
						outputfile << heu_actualUtility << " kept-alone solution " << evaluateUtilityBySol(type, depoCleanProb_keptalone, solKeptalone) << endl;

						resultfile << cplex_keptalone.getObjValue() << ","  << heu_actualUtility << "," << evaluateUtilityBySol(type, depoCleanProb_keptalone, solKeptalone) << ",";


						for (int i = 0; i < num_depo; i++)
							for (int j = 0; j < num_pool; j++)
								for (int k = 0; k < num_type; k++)
									solKeptalone[i][j][k] = 0;
						/*solKeptalone[3][0][0] = 19;
						solKeptalone[6][0][0] = 6;
						solKeptalone[1][1][1] = 25;
						solKeptalone[6][1][1] = 95;
						solKeptalone[5][2][0] = 83;
						solKeptalone[1][3][1] = 11;
						solKeptalone[4][3][1] = 38;
						solKeptalone[3][4][0] = 14;
						solKeptalone[1][5][2] = 58;
						solKeptalone[2][5][2] = 30;
						solKeptalone[0][6][1] = 73;
						solKeptalone[2][6][1] = 42;

						solKeptalone[1][0][1] = 26;
						solKeptalone[6][0][1] = 101;
						solKeptalone[5][1][0] = 83;
						solKeptalone[4][2][1] = 38;
						solKeptalone[2][2][1] = 17;
						solKeptalone[3][3][0] = 33;
						solKeptalone[1][4][2] = 50;
						solKeptalone[2][4][2] = 40;
						solKeptalone[0][5][1] = 25;
						solKeptalone[2][5][1] = 15;
						solKeptalone[0][6][1] = 48;
						solKeptalone[1][6][1] = 18;*/

						//env.out() << heu_actualUtility << " test " << evaluateUtilityBySol(type, depoCleanProb_keptalone, solKeptalone) << endl;

						//system("pause");
						for (i = 0; i < num_depo; i++) {
							for (j = 0; j < num_pool; j++) {
								delete[] solKeptalone[i][j];
							}
							delete[] solKeptalone[i];
						}
						delete[] solKeptalone;
						depoCleanProb_keptalone.clear();

						

						cplex_keptalone.end();
					}

					//Constraint: bound expected utility for initial feasiblity
					IloExpr v(env);
					v += expectedUtilityType[num_type];
					for (i = 0; i < num_depo; i++)
						for (j = 0; j < num_pool; j++)
							for (n = 0; n < num_type; n++)
								v += -type[n].unitUtility_def*allocateVolVars[i][j][n];
					model.add(v <= 0);
					v.end();



					model.add(IloMaximize(env, expectedUtilityType[num_type]));
					/*		if (minimize_num_pool) {
								IloExpr obj(env);
								obj += expectedUtilityType[num_pool];
								for (j = 1; j < num_pool; j++)
								obj -= cost_pool * IloSum(pooltypeBinVars[j]);
								model.add(IloMaximize(env, obj));
								obj.end();
								}
								else
								*/


					//env.out() << "reach here: after defining constraints" << endl;
					//		model.exportModel((filefolder + "model.lp").c_str());

					IloModel modelWithStart(env);
					modelWithStart.add(model);


					cplex.extract(model);
					cplex.exportModel((filefolder + "model.lp").c_str());

					IloCplex modelWithStartCplex(env);
					modelWithStartCplex.extract(modelWithStart);
					modelWithStartCplex.exportModel((filefolder + "model_start.lp").c_str());

					
					IloNumVarArray startVar(env);
					IloNumArray startVal(env);

					for (i = 0; i < num_depo; i++)
						for (j = 0; j < num_pool; j++){
							startVar.add(allocateVolVars[i][j]);
							startVal.add(allocateVolVarsSol_keptalone[i][j]);
							
						}
					modelWithStartCplex.addMIPStart(startVar, startVal);
					//system("pause");

					//ofstream LogFile(filefolder + "log.txt", ios_base::app);

					cplex.setOut(LogFile);//setout
					//cplex.setOut(env.getNullStream());
					//cplex.setWarning(env.getNullStream());


					cplex.setParam(IloCplex::TiLim, 600);
					cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.01);
					//cplex.setParam(IloCplex::VarSel,3);


					//env.out() << "reach here: before call  back" << endl;
					//cplex.use(LShapedLazyCallback(env, allocateVolVars, allocateBinVars, expectedUtilityType, max_pool_size, depoProb, typeUtility, cplex));
					cplex.use(LShapedLazyCallback(env, allocateVolVars, allocateBinVars, expectedUtilityType, depo, type, pools, cplex));
//					if (separateFracSols)
//						cplex.use(LShapedUserCallback(env, allocateVolVars, allocateBinVars, expectedUtilityType, depo, type, pools, cplex));
					//env.out() << "reach here: after call  back" << endl;

					IloNumArray2    ordpri_u(env, num_depo);
					IloNumArray2	ordpri_x(env, num_pool);
					if (set_priority_u) {
						for (i = 0; i < num_depo; i++) {
							ordpri_u[i] = IloNumArray(env, num_pool);
							for (j = 0; j < num_pool; j++)
								//ordpri_u[i][j] = (i *num_pool + j + 1) * 10;
								//ordpri_u[i][j] = (i + (num_pool - j)*num_depo) + set_priority_x*num_pool*num_type*2;
								ordpri_u[i][j] = (i*num_pool + (num_pool - j)) + set_priority_x * num_pool*num_type * 2;
							cplex.setPriorities(allocateBinVars[i], ordpri_u[i]);
						}
						//env.out() << "Set branch prioiry" << endl << endl;
					}

					if (set_priority_x) {
						for (j = 0; j < num_pool; j++) {
							ordpri_x[j] = IloNumArray(env, num_type);
							for (n = 0; n < num_type; n++)
								ordpri_x[j][n] = ((num_pool - j) * num_type + n + 1);
							cplex.setPriorities(pooltypeBinVars[j], ordpri_x[j]);
						}
					}

					cplex.solve();

					env.out() << "Solution status = " << cplex.getStatus() << endl;
					env.out() << "CPLEX status =    " << cplex.getCplexStatus() << endl;


					IloInt numNewMomIsolated = 0;
					IloInt numNewMomIsolatedStart = 0;


					if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible){

						env.out() << "The expected utility is " << cplex.getObjValue() << " " << cplex.getBestObjValue() << endl;
						env.out() << "The best utility is " << actualUtility_best << endl;
						outputfile << "The expected utility is " << cplex.getObjValue() << " " << cplex.getBestObjValue() << endl;
						outputfile << "The best utility is " << actualUtility_best << endl;
						resultfile << cplex.getObjValue() << "," << cplex.getBestObjValue() << "," << actualUtility_best;

						for (j = 0; j < num_pool; j++) {
							env.out() << "Pool " << j + 1 << " has ";
							outputfile << "Pool " << j + 1 << " has ";
							IloInt numDepoContained = 0;
							IloBool nonSplit = IloTrue;
							for (i = 0; i < num_depo; i++) {
								for (n = 0; n < num_type; n++) {
									if (cplex.getValue(allocateVolVars[i][j][n]) > RC_EPS) {
										env.out() << cplex.getValue(allocateVolVars[i][j][n]) << " from deposit " << i + 1 << " producing type " << n + 1 << endl;
										outputfile << cplex.getValue(allocateVolVars[i][j][n]) << " from deposit " << i + 1 << " producing type " << n + 1 << endl;
									}
									if (bestSolVol[i][j][n] > 1-RC_EPS) {
										numDepoContained++;
										if (bestSolVol[i][j][n]< depo[i].volume - RC_EPS || (!depo[i].new_mom))
											nonSplit = IloFalse;
									}
								}
							}
							if (nonSplit&&numDepoContained == 1) numNewMomIsolated++;
						}
						env.out() << "Number of isoloated new moms in risky pooling " << numNewMomIsolated << endl;
						outputfile << "Number of isoloated new moms in risky pooling " << numNewMomIsolated << endl;
						resultfile << numNewMomIsolated;
					}
					else{
						outputfile << "The expected utility is (cannot find a feasible soln) " << " " << cplex.getBestObjValue() << endl;
						outputfile << "The best utility is  (cannot find a feasible soln)  " << endl;
						resultfile << "noSoln," << "noSoln,NA" ;
					
					}

					for (j = 0; j < num_pool; j++) {
						env.out() << "Pool " << j + 1 << " has ";
						outputfile << "Pool " << j + 1 << " has ";
						for (i = 0; i < num_depo; i++) {
							for (n = 0; n < num_type; n++) {
								if (bestSolVol[i][j][n] > RC_EPS) {
									env.out() << bestSolVol[i][j][n] << " from deposit " << i + 1 << " producing type " << n + 1 << endl;
									outputfile << bestSolVol[i][j][n] << " from deposit " << i + 1 << " producing type " << n + 1 << endl;
								}
							}
						}
					}


					modelWithStartCplex.setOut(LogFile);//setout
					//cplex.setOut(env.getNullStream());
					//cplex.setWarning(env.getNullStream());


					modelWithStartCplex.setParam(IloCplex::TiLim, 600);
					modelWithStartCplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.01);

					modelWithStartCplex.use(LShapedLazyCallback(env, allocateVolVars, allocateBinVars, expectedUtilityType, depo, type, pools, modelWithStartCplex));

					modelWithStartCplex.solve();
					if (set_priority_u) {
						for (i = 0; i < num_depo; i++) {
							ordpri_u[i] = IloNumArray(env, num_pool);
							for (j = 0; j < num_pool; j++)
								ordpri_u[i][j] = (i*num_pool + (num_pool - j)) + set_priority_x * num_pool*num_type * 2;
							modelWithStartCplex.setPriorities(allocateBinVars[i], ordpri_u[i]);
						}
					}

					if (set_priority_x) {
						for (j = 0; j < num_pool; j++) {
							ordpri_x[j] = IloNumArray(env, num_type);
							for (n = 0; n < num_type; n++)
								ordpri_x[j][n] = ((num_pool - j) * num_type + n + 1);
							modelWithStartCplex.setPriorities(pooltypeBinVars[j], ordpri_x[j]);
						}
					}

					ordpri_u.end();
					ordpri_x.end();

					env.out() << "The expected utility with a start " << modelWithStartCplex.getObjValue() << " " << modelWithStartCplex.getBestObjValue() << endl;
					outputfile << "The expected utility with a start " << modelWithStartCplex.getObjValue() << " " << modelWithStartCplex.getBestObjValue() << endl;

					resultfile << ",," << modelWithStartCplex.getObjValue() << " ," << modelWithStartCplex.getBestObjValue() << " ,";

					for (j = 0; j < num_pool; j++) {
						env.out() << "Pool " << j + 1 << " has ";
						outputfile << "Pool " << j + 1 << " has ";
						IloInt numDepoContainedStart = 0;
						IloBool nonSplitStart = IloTrue;
						for (i = 0; i < num_depo; i++) {
							for (n = 0; n < num_type; n++) {
								if (modelWithStartCplex.getValue(allocateVolVars[i][j][n]) > RC_EPS ) {
									numDepoContainedStart++;
									if (modelWithStartCplex.getValue(allocateVolVars[i][j][n]) < depo[i].volume - RC_EPS || (!depo[i].new_mom)){
										nonSplitStart = IloFalse;
										//env.out() << i + 1 << " deposit " << (modelWithStartCplex.getValue(allocateVolVars[i][j][n]) < depo[i].volume - RC_EPS) << " " << (!depo[i].new_mom) << endl;
									}
									env.out() << modelWithStartCplex.getValue(allocateVolVars[i][j][n]) << " (with a start) from deposit " << i + 1 << " producing type " << n + 1 << endl;
									outputfile << modelWithStartCplex.getValue(allocateVolVars[i][j][n]) << " (with a start) from deposit " << i + 1 << " producing type " << n + 1 << endl;
								}

							}
						}						
						if (nonSplitStart&&numDepoContainedStart == 1) numNewMomIsolatedStart++;
					}

					env.out() << "Number of isoloated new moms with a start " << numNewMomIsolatedStart << endl;
					outputfile << "Number of isoloated new moms with a start  " << numNewMomIsolatedStart << endl;

					env.out() << "Number of cuts added " << cut_count << endl;
					outputfile << "Number of cuts added " << cut_count << endl;

					resultfile << numNewMomIsolatedStart << ",";

					vector<double> depoCleanProb_double;
					for (i = 0; i < num_depo; i++)
						depoCleanProb_double.push_back(depo[i].bacteriaFreeProb);
					double utilityReevaluated = evaluateUtilityBySol(type, depoCleanProb_double, bestSolVol);
					env.out() << "separate evaluate best solution: " << utilityReevaluated << endl;
					outputfile<< "separate evaluate best solution: " << utilityReevaluated << endl;


					double utilityOtherSetting = evaluateUtilityBySol(type, depoCleanProb_double, bestSolVolOtherSetting);
					env.out() << "other best solution: " << utilityOtherSetting << endl;
					outputfile << "other best solution: " << utilityOtherSetting << endl;

					

					IloInt numNewMomIsolatedOther = 0;
					if (utilityOtherSetting < utilityReevaluated) {
						for (i = 0; i < num_depo; i++)
							for (j = 0; j < num_pool; j++)
								for (n = 0; n < num_type; n++){
									bestSolVolOtherSetting[i][j][n] = bestSolVol[i][j][n];
								}
						resultfile  << "NA ,";
					}
					else {
						for (j = 0; j < num_pool; j++) {
							env.out() << "Pool " << j + 1 << " has ";
							outputfile << "Pool " << j + 1 << " has ";
							IloInt numDepoContainedOther = 0;
							IloInt nonSplitOther = IloTrue;
							for (i = 0; i < num_depo; i++) {
								for (n = 0; n < num_type; n++) {
									if (bestSolVolOtherSetting[i][j][n] > RC_EPS) {
										env.out() << bestSolVolOtherSetting[i][j][n] << " from deposit " << i + 1 << " producing type " << n + 1 << endl;
										outputfile << bestSolVolOtherSetting[i][j][n] << " from deposit " << i + 1 << " producing type " << n + 1 << endl;
										numDepoContainedOther++;
										if (bestSolVolOtherSetting[i][j][n] < depo[i].volume - RC_EPS || (!depo[i].new_mom)){
											nonSplitOther = IloFalse;
										}										
									}
								}
							}
						
						if (nonSplitOther&&numDepoContainedOther == 1) numNewMomIsolatedOther++;
						//env.out() << numNewMomIsolatedOther << " isolated other setting" << endl;
						
						}
						resultfile << numNewMomIsolatedOther << ",";
					}



					resultfile << utilityOtherSetting << "," << endl;

					for (i = 0; i < num_depo; i++) {
						for (j = 0; j < num_pool; j++){
							//for (n = 0; n < num_type; n++)
							//if (allocateVolVarsSol[i][j][n] > RC_EPS)
							//env.out() << "Deposit " << i + 1 << " Pool " << j + 1 << " producing " << n + 1 << " Volume is " << allocateVolVarsSol[i][j][n] << endl;
							allocateVolVarsSol_keptalone[i][j].end();
						}
						allocateVolVarsSol_keptalone[i].end();
					}
					allocateVolVarsSol_keptalone.end();

					//LogFile.close();
					//system("pause");

					for (i = 0; i < num_depo; i++) {
						for (j = 0; j < num_pool; j++) {
							delete[] bestSolVol[i][j];
						}
						delete[] bestSolVol[i];
					}
					delete [] bestSolVol;
									if (m == num_established && mm == num_newmom) {
					for (int i = 0; i < num_depo; i++) {
						for (int j = 0; j < num_pool; j++) {
							delete[] bestSolVolOtherSetting[i][j];
						}
						delete[] bestSolVolOtherSetting[i];
					}
					delete[] bestSolVolOtherSetting;
				}
				}

				catch (IloException& e) {
					cerr << "Concert exception caught: " << e << endl;
					//system("pause");
				}

				env.end();

				depo.clear();
				type.clear();
				for (int j = 0; j < pools.size(); j++)
					pools[j].depo_idx.clear();
				pools.clear();


			}
		}

	}
LogFile.close();
outputfile.close();
resultfile.close();
system("pause");
	return 0;

}  // END main

//IloBool OptCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const IloNumArray depoProb, const IloNumArray typeUtility,  IloCplex cplex, IloExpr v_singlecut, IloExpr2dArray v_multicut, IloNum& singlecutRHS, IloNumArray multicutRHS)
IloBool OptSingleCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const vector<deposit> depo, const vector<milktype> type, IloCplex cplex, IloExpr v_singlecut, IloNum& singlecutRHS, vector<pool> pools)
//IloBool OptCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const vector<deposit> depo, const vector<milktype> type, IloCplex cplex, IloArray<IloExpr> v_multicut, IloNumArray singlecutRHS, vector<pool> pools)
{
	IloBool violatedCutFound = IloFalse;

	IloEnv env = cplex.getEnv();
	IloModel mod = cplex.getModel();

	IloInt num_depo = allocateVolVars.getSize();
	IloInt num_pool = allocateVolVars[0].getSize();
	IloInt num_type = allocateVolVars[0][0].getSize();
	IloInt i, i1, j, j1, n;

	
	//IloNum sumUnitUtility = IloSum(typeUtility);

		IloNum actualUtility = 0.0;
		v_singlecut.clear();
		singlecutRHS = 0.0;

		ObtainPoolInfo(allocateVolVarsSol, depo, pools);
		IloNumArray poolCleanProb(env, num_pool);
		for (j = 0; j < num_pool; j++)
			poolCleanProb[j] = pools[j].clean_prob;
		IloNumArray productionTargets(env, num_type);
		IloNumArray depoCleanProb(env, num_depo);
		IloNumArray depoVol(env, num_depo);

		//for (j = 0; j < num_pool; j++){
			//env.out() << j + 1 << "contains depos ";
			//for (i = 0; i < pools[j].depo_idx.size(); i++)
			//	env.out() << pools[j].depo_idx[i] + 1 << " ";
			//env.out() << " produce type " << pools[j].type_idx + 1 << " with prob " << pools[j].clean_prob << endl;
		//}

		//IloNumArray2 TypePoolVol(env, num_type);
		//IloNumArray2 TypePoolProb(env, num_type);
		//for (n = 0; n < num_type; n++){
		//	TypePoolVol[n] = IloNumArray(env);
		//	TypePoolProb[n] = IloNumArray(env);
		//}
		//IloBoolArray producedType = ObtainTypeProduction(TypePoolVol, TypePoolProb, pools);
		for (i = 0; i < num_depo; i++) {
			depoCleanProb[i] = depo[i].bacteriaFreeProb;
			depoVol[i] = depo[i].volume;
		}
		

		IloNum bigM = max_pool_size;
		IloNumArray2 bigMArray(env, num_pool);
		for (j = 0; j < num_pool; j++) {
			bigMArray[j] = IloNumArray(env, num_type);
			for (n = 0; n < num_type; n++) {
				bigMArray[j][n] = 0;
				for (i = 0; i < num_depo; i++)
					if(allocateVolVarsSol[i][j][n]>RC_EPS)
						bigMArray[j][n] += allocateVolVarsSol[i][j][n];
			}
		}
		for (int n = 0; n < num_type; n++)
			productionTargets[n] = type[n].productionTarget;

		IloNumArray meetTargetProbAll = CalculateMeetTargetProb(pools, productionTargets, depoCleanProb, depoVol);
		
		for (int n = 0; n < num_type; n++){
			singlecutRHS -= type[n].productionTarget * (type[n].unitUtility_def - type[n].unitUtility_sur)*meetTargetProbAll[n];
				//env.out() << "Type " << n + 1 << "Prob" << meetTargetProbAll[n] << " " << singlecutRHS << endl;

		}
		//env.out() << singlecutRHS+1;
		actualUtility -= singlecutRHS;
		//env.out() << " RHS: actual utility" << actualUtility << endl;

			

		for (j = 0; j < num_pool; j++) {
			IloBool emptyPool = 1;
			IloNum poolNotContaminationProb = 1;
			IloNum depositsHaventContaminate = 1;// The probability that the contained deposits haven't conataminated the pool yet.
			IloNum poolUtility = 0;
			IloNumArray depo_prob_temp(env, num_depo);
			for (i = 0; i < num_depo; i++){
				depo_prob_temp[i] = depo[i].bacteriaFreeProb;
				//env.out() << depo_prob_temp[i] << endl;
			}
			if (pools[j].is_used)
				for (i = 0; i < pools[j].depo_idx.size(); i++){
					depo_prob_temp[pools[j].depo_idx[i]] = 1;
				}
			//env.out() << depo_prob_temp[0] << " " << depo_prob_temp[1] << " " << depo_prob_temp[2] << endl;
			//for (n = 0; n < num_type; n++) env.out() << productionTarget[n] << endl;
			IloNumArray meetTargetConProbAll_Pool = CalculateMeetTargetProb(pools, productionTargets, depo_prob_temp, depoVol);
			for (int n = 0; n < num_type; n++){
				//env.out() << "meetTargetConProbAll_Pool " << meetTargetConProbAll_Pool << endl;
			for (i = 0; i < num_depo; i++){
				v_singlecut += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVars[i][j][n];
					if (allocateVolVarsSol[i][j][n]>RC_EPS){
						actualUtility += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVarsSol[i][j][n];
						//env.out() << "actual utility increment " << (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] << " " << allocateVolVarsSol[i][j][n] <<" " << (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVarsSol[i][j][n] << " "<< actualUtility<<endl;
						//env.out() << "actual utility increment ------ " << n <<" " << type[n].unitUtility_def<< " " << (1 - meetTargetConProbAll_Pool[n]) << " " << type[n].unitUtility_sur << " " << (meetTargetConProbAll_Pool[n]) << " " << poolCleanProb[j] << endl;
					//env.out() << i << j<< n << "actual utility" << actualUtility << endl;
					}
					}
			//env.out() << v_singlecut << endl;
			//if(meetTargetProbAll[n]<1)
			singlecutRHS -= bigMArray[j][n] * (type[n].unitUtility_def*(-(1-meetTargetConProbAll_Pool[n]) * poolCleanProb[j] +1- meetTargetProbAll[n]) + type[n].unitUtility_sur*(meetTargetProbAll[n] - meetTargetConProbAll_Pool[n] * poolCleanProb[j]));
			//else singlecutRHS -= bigM * type[n].unitUtility_sur;
//			env.out() << "Pool " << j <<" is clean implying" << n  << " : " << productionTarget[n] << "probability of meeting target " << meetTargetConProbAll_Pool[n] << " " << meetTargetProbAll[n] << endl;
			//env.out() << bigM * (type[n].unitUtility_def*(-(1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j] + 1 - meetTargetProbAll[n]) + type[n].unitUtility_sur*(meetTargetProbAll[n] - meetTargetConProbAll_Pool[n] * poolCleanProb[j])) << " singlecutRHS "<< singlecutRHS << endl;

			//env.out() << " singlecutRHS "<< (1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j]  << " " << 1 - meetTargetProbAll[n] << "------ " << (-(1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j] + 1 - meetTargetProbAll[n])  << "  "<< (meetTargetProbAll[n] - meetTargetConProbAll_Pool[n] * poolCleanProb[j]) << " singlecutRHS " <<  singlecutRHS << endl;
			}

			for (i = 0; i < num_depo; i++){
				depo_prob_temp[i] = depo[i].bacteriaFreeProb;
			}
			if (pools[j].is_used){
				for (i = 0; i < pools[j].depo_idx.size(); i++){
					IloNum isContaminatedDepo = 1;
					for (i1 = 0; i1 < i; i1++){
						isContaminatedDepo *= depo[pools[j].depo_idx[i1]].bacteriaFreeProb;
							depo_prob_temp[pools[j].depo_idx[i1]] = 1;
					}
					isContaminatedDepo *= 1-depo[pools[j].depo_idx[i]].bacteriaFreeProb;
					depo_prob_temp[pools[j].depo_idx[i]] = 0;
					IloNumArray meetTargetConProbAll_Depo = CalculateMeetTargetProb(pools, productionTargets, depo_prob_temp, depoVol);
					IloNum b_ij = 0;
						for (IloInt n1 = 0; n1 < num_type; n1++){
							b_ij += bigMArray[j][n1] *(type[n1].unitUtility_def*(1 - meetTargetConProbAll_Depo[n1]) + type[n1].unitUtility_sur*(meetTargetConProbAll_Depo[n1]))*isContaminatedDepo;
						}
						v_singlecut -= b_ij * allocateBinVars[pools[j].depo_idx[i]][j];
						meetTargetConProbAll_Depo.end();
				}
			}
			depo_prob_temp.end();
			meetTargetConProbAll_Pool.end();
			/*
			for (i = 0; i < num_depo; i++) {
				if (allocateBinVarsSol[i][j] >= 1 - RC_EPS) {
					emptyPool = 0;
					v_singlecut -= max_pool_size*sumUnitUtility*poolNotContaminationProb*depoProb[i]* allocateBinVars[i][j];
					//for (j1 = 0; j1<num_pool; j1++)
					//	v_multicut[j][j1] -= max_pool_size*sumUnitUtility*poolNotContaminationProb*depoProb[i]* allocateBinVars[i][j1];
					poolNotContaminationProb *= (1 - depoProb[i]);
				}
				for (n = 0; n < num_type; n++)
					if (allocateVolVarsSol[i][j][n] >= RC_EPS) {
						poolUtility += allocateVolVarsSol[i][j][n]*typeUtility[n];
					}
			}
			for (i = 0; i < num_depo; i++) {
				for (n = 0; n < num_type; n++) {
					v_singlecut += typeUtility[n] *poolNotContaminationProb * allocateVolVars[i][j][n];
					//for (j1 = 0; j1<num_pool; j1++)
					//	v_multicut[j][j1] += typeUtility[n] *poolNotContaminationProb * allocateVolVars[i][j1][n];
				}
			}*/
			//actualUtility += poolUtility*poolNotContaminationProb;
			//for (n = 0; n < num_type; n++) 
			//	optCutRHS -= max_pool_size*(1 - poolNotContaminationProb)*type[n].unitUtility;
			//singlecutRHS -= max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
			//multicutRHS[j] = -max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
			//for (j1 = 0; j1<num_pool; j1++)
			//	v_multicut[j][j1] += max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
						//env.out() << v_multicut[j][0] << " " << multicutRHS[j] << endl;
		}

			/*for (j = 0; j<num_pool; j++)
				for (j1 = 0; j1<num_pool; j1++) {
					v_multicut[j][j1] -= expectedUtilityType[j1];
				}
				*/
			v_singlecut -= expectedUtilityType[num_type];
			//env.out() << "actual utility" << actualUtility << ". best utility" << actualUtility_best << endl;
			if (actualUtility_best < actualUtility) {
				actualUtility_best = actualUtility;
				for (i = 0; i < num_depo; i++) {
					for (j = 0; j < num_pool; j++) {
						for (n = 0; n < num_type; n++) {
							bestSolVol[i][j][n] = allocateVolVarsSol[i][j][n];
							//cout << bestSolVol[i][j][n];
							//system("pause");
						}
					}
				}
			}
			if (actualUtility < expectedUtilityTypeSol[num_type] - RC_EPS){
				violatedCutFound = IloTrue;
			//system("pause");
				cut_count++;
				if (cut_count % 100 == 0) 
					env.out() << "Have added cuts ---- " << cut_count << endl;
			}
			for (j = 0; j < num_pool; j++)
				bigMArray[j].end();
			bigMArray.end();
			depoCleanProb.end();
			poolCleanProb.end();



			//env.out() << actualUtility_best  << endl;
			//poolCleanProb.end();
			//productionTargets.end();
			//depoCleanProb.end();
	return violatedCutFound;

} // END OptCutAdded

IloBool OptAdditionalCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const vector<deposit> depo, const vector<milktype> type, IloCplex cplex, IloArray<IloExpr> v_multicut, IloNum& multicutRHS, vector<pool> pools, IloBool & bestUtilityImproved)
{
	IloBool violatedCutFound = IloFalse;

	IloEnv env = cplex.getEnv();
	IloModel mod = cplex.getModel();

	IloInt num_depo = allocateVolVars.getSize();
	IloInt num_pool = allocateVolVars[0].getSize();
	IloInt num_type = allocateVolVars[0][0].getSize();
	IloInt i, i1, j, j1, n;

	//IloNum sumUnitUtility = IloSum(typeUtility);

	IloNum actualUtility = 0.0;
	//v_singlecut.clear();
	multicutRHS = 0.0;

	ObtainPoolInfo(allocateVolVarsSol, depo, pools);
	IloNumArray poolCleanProb(env, num_pool);
	for (j = 0; j < num_pool; j++)
		poolCleanProb[j] = pools[j].clean_prob;
	IloNumArray productionTargets(env, num_type);
	IloNumArray depoCleanProb(env, num_depo);
	IloNumArray depoVol(env, num_depo);

	//for (j = 0; j < num_pool; j++){
		//env.out() << j + 1 << "contains depos ";
		//for (i = 0; i < pools[j].depo_idx.size(); i++)
		//	env.out() << pools[j].depo_idx[i] + 1 << " ";
		//env.out() << " produce type " << pools[j].type_idx + 1 << " with prob " << pools[j].clean_prob << endl;
	//}

	//IloNumArray2 TypePoolVol(env, num_type);
	//IloNumArray2 TypePoolProb(env, num_type);
	//for (n = 0; n < num_type; n++){
	//	TypePoolVol[n] = IloNumArray(env);
	//	TypePoolProb[n] = IloNumArray(env);
	//}
	//IloBoolArray producedType = ObtainTypeProduction(TypePoolVol, TypePoolProb, pools);
	for (i = 0; i < num_depo; i++) {
		depoCleanProb[i] = depo[i].bacteriaFreeProb;
		depoVol[i] = depo[i].volume;
	}


	IloNum bigM = max_pool_size;
	IloNumArray2 bigMArray(env, num_pool);
	for (j = 0; j < num_pool; j++) {
		bigMArray[j] = IloNumArray(env, num_type);
		for (n = 0; n < num_type; n++) {
			bigMArray[j][n] = 0;
			for (i = 0; i < num_depo; i++)
				if (allocateVolVarsSol[i][j][n] > RC_EPS)
					bigMArray[j][n] += allocateVolVarsSol[i][j][n];
		}
	}
	for (int n = 0; n < num_type; n++)
		productionTargets[n] = type[n].productionTarget;

	IloNumArray meetTargetProbAll = CalculateMeetTargetProb(pools, productionTargets, depoCleanProb, depoVol);

	for (int n = 0; n < num_type; n++) {
		multicutRHS -= type[n].productionTarget * (type[n].unitUtility_def - type[n].unitUtility_sur)*meetTargetProbAll[n];
		//env.out() << "Type " << n + 1 << "Prob" << meetTargetProbAll[n] << " " << singlecutRHS << endl;

	}
	//env.out() << singlecutRHS+1;
	actualUtility -= multicutRHS;
	//env.out() << " RHS: actual utility" << actualUtility << endl;



	for (j = 0; j < num_pool; j++) {
		IloBool emptyPool = 1;
		IloNum poolNotContaminationProb = 1;
		IloNum depositsHaventContaminate = 1;// The probability that the contained deposits haven't conataminated the pool yet.
		IloNum poolUtility = 0;
		IloNumArray depo_prob_temp(env, num_depo);
		for (i = 0; i < num_depo; i++) {
			depo_prob_temp[i] = depo[i].bacteriaFreeProb;
			//env.out() << depo_prob_temp[i] << endl;
		}
		if (pools[j].is_used)
			for (i = 0; i < pools[j].depo_idx.size(); i++) {
				depo_prob_temp[pools[j].depo_idx[i]] = 1;
			}
		IloNumArray meetTargetConProbAll_Pool = CalculateMeetTargetProb(pools, productionTargets, depo_prob_temp, depoVol);
		for (int n = 0; n < num_type; n++) {
			for (i = 0; i < num_depo; i++) {
				v_multicut[0] += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVars[i][j][n];
				for (j1 = 1; j1 < num_pool; j1++) {
					if(j1==j) v_multicut[j1] += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j-1] * allocateVolVars[i][j-1][n];
					else if (j1 == j+1) v_multicut[j1] += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j + 1] * allocateVolVars[i][j + 1][n];
					else v_multicut[j1] += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVars[i][j][n];
				}

				if (allocateVolVarsSol[i][j][n] > RC_EPS) {
					actualUtility += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVarsSol[i][j][n];
				}
			}
			multicutRHS -= bigMArray[j][n] * (type[n].unitUtility_def*(-(1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j] + 1 - meetTargetProbAll[n]) + type[n].unitUtility_sur*(meetTargetProbAll[n] - meetTargetConProbAll_Pool[n] * poolCleanProb[j]));
		}

		for (i = 0; i < num_depo; i++) {
			depo_prob_temp[i] = depo[i].bacteriaFreeProb;
		}
		if (pools[j].is_used) {
			for (i = 0; i < pools[j].depo_idx.size(); i++) {
				IloNum isContaminatedDepo = 1;
				for (i1 = 0; i1 < i; i1++) {
					isContaminatedDepo *= depo[pools[j].depo_idx[i1]].bacteriaFreeProb;
					depo_prob_temp[pools[j].depo_idx[i1]] = 1;
				}
				isContaminatedDepo *= 1 - depo[pools[j].depo_idx[i]].bacteriaFreeProb;
				depo_prob_temp[pools[j].depo_idx[i]] = 0;
				IloNumArray meetTargetConProbAll_Depo = CalculateMeetTargetProb(pools, productionTargets, depo_prob_temp, depoVol);
				IloNum b_ij = 0;

				for (IloInt n1 = 0; n1 < num_type; n1++) {
					b_ij += bigMArray[j][n1] * (type[n1].unitUtility_def*(1 - meetTargetConProbAll_Depo[n1]) + type[n1].unitUtility_sur*(meetTargetConProbAll_Depo[n1]))*isContaminatedDepo;
				}
				v_multicut[0] -= b_ij * allocateBinVars[pools[j].depo_idx[i]][j];
				for (j1 = 1; j1 < num_pool; j1++) {
					if (j1 == j) v_multicut[j1] -= b_ij * allocateBinVars[pools[j].depo_idx[i]][j-1];
					else if (j1 == j + 1) v_multicut[j1] -= b_ij * allocateBinVars[pools[j].depo_idx[i]][j+1];
					else v_multicut[j1] -= b_ij * allocateBinVars[pools[j].depo_idx[i]][j];
				}

			}
		}
		depo_prob_temp.end();
		/*
		for (i = 0; i < num_depo; i++) {
			if (allocateBinVarsSol[i][j] >= 1 - RC_EPS) {
				emptyPool = 0;
				v_singlecut -= max_pool_size*sumUnitUtility*poolNotContaminationProb*depoProb[i]* allocateBinVars[i][j];
				//for (j1 = 0; j1<num_pool; j1++)
				//	v_multicut[j][j1] -= max_pool_size*sumUnitUtility*poolNotContaminationProb*depoProb[i]* allocateBinVars[i][j1];
				poolNotContaminationProb *= (1 - depoProb[i]);
			}
			for (n = 0; n < num_type; n++)
				if (allocateVolVarsSol[i][j][n] >= RC_EPS) {
					poolUtility += allocateVolVarsSol[i][j][n]*typeUtility[n];
				}
		}
		for (i = 0; i < num_depo; i++) {
			for (n = 0; n < num_type; n++) {
				v_singlecut += typeUtility[n] *poolNotContaminationProb * allocateVolVars[i][j][n];
				//for (j1 = 0; j1<num_pool; j1++)
				//	v_multicut[j][j1] += typeUtility[n] *poolNotContaminationProb * allocateVolVars[i][j1][n];
			}
		}*/
		//actualUtility += poolUtility*poolNotContaminationProb;
		//for (n = 0; n < num_type; n++) 
		//	optCutRHS -= max_pool_size*(1 - poolNotContaminationProb)*type[n].unitUtility;
		//singlecutRHS -= max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
		//multicutRHS[j] = -max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
		//for (j1 = 0; j1<num_pool; j1++)
		//	v_multicut[j][j1] += max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
					//env.out() << v_multicut[j][0] << " " << multicutRHS[j] << endl;
	}

	/*for (j = 0; j<num_pool; j++)
		for (j1 = 0; j1<num_pool; j1++) {
			v_multicut[j][j1] -= expectedUtilityType[j1];
		}
		*/
	for(j = 0; j< num_pool; j++)
		v_multicut[j] -= expectedUtilityType[num_type];
	//env.out() << "actual utility" << actualUtility << ". expected utility" << expectedUtilityTypeSol[num_pool] << endl;
	if (actualUtility_best < actualUtility) {
		bestUtilityImproved = IloTrue;
		actualUtility_best = actualUtility;
		for (j = 0; j < num_pool; j++) {
			for (i = 0; i < num_depo; i++) {
				for (n = 0; n < num_type; n++) {
					bestSolVol[i][j][n] = allocateVolVarsSol[i][j][n];
					//if (allocateVolVarsSol[i][j][n] > RC_EPS)
						//env.out() << "Deposit " << i + 1 << " to Pool " << j + 1 << " producing type " << n + 1 << " Volume is " << allocateVolVarsSol[i][j][n] << endl;
				}
			}
		}

	}
	if (actualUtility < expectedUtilityTypeSol[num_type] - RC_EPS) {
		violatedCutFound = IloTrue;
		//system("pause");

		//env.out() << v_multicut[0] << " >= " << multicutRHS << endl;
		//env.out() << evaluateActualUtility(depoCleanProb, pools, type) << endl;
		//system("pause");
		
		cut_count++;
		if (cut_count % 100 == 0){
			env.out() << "Have added cuts ---- " << cut_count << endl;
			//env.out() << "best actual utility" << actualUtility_best << ". expected utility" << expectedUtilityTypeSol[num_type] << endl;
		}
	}
	for (j = 0; j < num_pool; j++)
		bigMArray[j].end();
	bigMArray.end();
	depoCleanProb.end();
	poolCleanProb.end();
	depoVol.end();
	productionTargets.end();



	//env.out() << actualUtility_best  << endl;
	//poolCleanProb.end();
	//productionTargets.end();
	//depoCleanProb.end();
	return violatedCutFound;
} // END OptCutAdded


IloBool OptMultiCutAdded(const IloNumVar3dArray allocateVolVars, const IloNumArray3 allocateVolVarsSol, const IloNumVar2dArray allocateBinVars, const IloNumArray2 allocateBinVarsSol, const IloNumVarArray expectedUtilityType, const IloNumArray expectedUtilityTypeSol, const IloNum max_pool_size, const vector<deposit> depo, const vector<milktype> type, IloCplex cplex, IloArray<IloExpr> v_multicut, IloNumArray multicutRHS, vector<pool> pools)
{
	IloBool violatedCutFound = IloFalse;

	IloEnv env = cplex.getEnv();
	IloModel mod = cplex.getModel();

	IloInt num_depo = allocateVolVars.getSize();
	IloInt num_pool = allocateVolVars[0].getSize();
	IloInt num_type = allocateVolVars[0][0].getSize();
	IloInt i, i1, j, j1, n;


	//IloNum sumUnitUtility = IloSum(typeUtility);

	IloNum actualUtility = 0.0;
	for (n = 0; n < num_type; n++) {
		v_multicut[n].clear();
		multicutRHS[n] = 0.0;
	}

	ObtainPoolInfo(allocateVolVarsSol, depo, pools);
	IloNumArray poolCleanProb(env, num_pool);
	for (j = 0; j < num_pool; j++)
		poolCleanProb[j] = pools[j].clean_prob;
	IloNumArray productionTargets(env, num_type);
	IloNumArray depoCleanProb(env, num_depo);
	IloNumArray depoVol(env, num_depo);
	for (i = 0; i < num_depo; i++) {
		depoCleanProb[i] = depo[i].bacteriaFreeProb;
		depoVol[i] = depo[i].volume;
	}

	IloNum bigM = max_pool_size;
	IloNumArray2 bigMArray(env, num_pool);
	for (j = 0; j < num_pool; j++) {
		bigMArray[j] = IloNumArray(env, num_type);
		for (n = 0; n < num_type; n++) {
			bigMArray[j][n] = 0;
			for (i = 0; i < num_depo; i++)
				if (allocateVolVarsSol[i][j][n] > RC_EPS)
					bigMArray[j][n] += allocateVolVarsSol[i][j][n];
		}
	}
	for (int n = 0; n < num_type; n++)
		productionTargets[n] = type[n].productionTarget;

	IloNumArray meetTargetProbAll = CalculateMeetTargetProb(pools, productionTargets, depoCleanProb, depoVol);

	for (int n = 0; n < num_type; n++) {
		multicutRHS[n] -= type[n].productionTarget * (type[n].unitUtility_def - type[n].unitUtility_sur)*meetTargetProbAll[n];
		//env.out() << "Type " << n + 1 << "Prob" << meetTargetProbAll[n] << " " << singlecutRHS << endl;
		actualUtility -= multicutRHS[n];
	}
	//env.out() << singlecutRHS+1;

	//env.out() << " RHS: actual utility" << actualUtility << endl;



	for (j = 0; j < num_pool; j++) {
		IloBool emptyPool = 1;
		IloNum poolNotContaminationProb = 1;
		IloNum depositsHaventContaminate = 1;// The probability that the contained deposits haven't conataminated the pool yet.
		IloNum poolUtility = 0;
		IloNumArray depo_prob_temp(env, num_depo);
		for (i = 0; i < num_depo; i++) {
			depo_prob_temp[i] = depo[i].bacteriaFreeProb;
			//env.out() << depo_prob_temp[i] << endl;
		}
		if (pools[j].is_used)
			for (i = 0; i < pools[j].depo_idx.size(); i++) {
				depo_prob_temp[pools[j].depo_idx[i]] = 1;
			}
		//env.out() << depo_prob_temp[0] << " " << depo_prob_temp[1] << " " << depo_prob_temp[2] << endl;
		//for (n = 0; n < num_type; n++) env.out() << productionTarget[n] << endl;
		IloNumArray meetTargetConProbAll_Pool = CalculateMeetTargetProb(pools, productionTargets, depo_prob_temp, depoVol);
		for (int n = 0; n < num_type; n++) {
			//env.out() << "meetTargetConProbAll_Pool " << meetTargetConProbAll_Pool << endl;
			for (i = 0; i < num_depo; i++) {
				v_multicut[n] += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVars[i][j][n];
				if (allocateVolVarsSol[i][j][n] > RC_EPS) {
					actualUtility += (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVarsSol[i][j][n];
					//env.out() << "actual utility increment " << (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] << " " << allocateVolVarsSol[i][j][n] <<" " << (type[n].unitUtility_def*(1 - meetTargetConProbAll_Pool[n]) + type[n].unitUtility_sur*(meetTargetConProbAll_Pool[n]))*poolCleanProb[j] * allocateVolVarsSol[i][j][n] << " "<< actualUtility<<endl;
					//env.out() << "actual utility increment ------ " << n <<" " << type[n].unitUtility_def<< " " << (1 - meetTargetConProbAll_Pool[n]) << " " << type[n].unitUtility_sur << " " << (meetTargetConProbAll_Pool[n]) << " " << poolCleanProb[j] << endl;
				//env.out() << i << j<< n << "actual utility" << actualUtility << endl;
				}
			}
			//env.out() << v_singlecut << endl;
			//if(meetTargetProbAll[n]<1)
			multicutRHS[n] -= bigMArray[j][n] * (type[n].unitUtility_def*(-(1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j] + 1 - meetTargetProbAll[n]) + type[n].unitUtility_sur*(meetTargetProbAll[n] - meetTargetConProbAll_Pool[n] * poolCleanProb[j]));
			//else singlecutRHS -= bigM * type[n].unitUtility_sur;
//			env.out() << "Pool " << j <<" is clean implying" << n  << " : " << productionTarget[n] << "probability of meeting target " << meetTargetConProbAll_Pool[n] << " " << meetTargetProbAll[n] << endl;
			//env.out() << bigM * (type[n].unitUtility_def*(-(1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j] + 1 - meetTargetProbAll[n]) + type[n].unitUtility_sur*(meetTargetProbAll[n] - meetTargetConProbAll_Pool[n] * poolCleanProb[j])) << " singlecutRHS "<< singlecutRHS << endl;

			//env.out() << " singlecutRHS "<< (1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j]  << " " << 1 - meetTargetProbAll[n] << "------ " << (-(1 - meetTargetConProbAll_Pool[n]) * poolCleanProb[j] + 1 - meetTargetProbAll[n])  << "  "<< (meetTargetProbAll[n] - meetTargetConProbAll_Pool[n] * poolCleanProb[j]) << " singlecutRHS " <<  singlecutRHS << endl;
		}

		for (i = 0; i < num_depo; i++) {
			depo_prob_temp[i] = depo[i].bacteriaFreeProb;
		}
		if (pools[j].is_used) {
			for (i = 0; i < pools[j].depo_idx.size(); i++) {
				IloNum isContaminatedDepo = 1;
				for (i1 = 0; i1 < i; i1++) {
					isContaminatedDepo *= depo[pools[j].depo_idx[i1]].bacteriaFreeProb;
					depo_prob_temp[pools[j].depo_idx[i1]] = 1;
				}
				isContaminatedDepo *= 1 - depo[pools[j].depo_idx[i]].bacteriaFreeProb;
				depo_prob_temp[pools[j].depo_idx[i]] = 0;
				IloNumArray meetTargetConProbAll_Depo = CalculateMeetTargetProb(pools, productionTargets, depo_prob_temp, depoVol);

				for (IloInt n1 = 0; n1 < num_type; n1++) {
					IloNum b_ij = bigMArray[j][n1] * (type[n1].unitUtility_def*(1 - meetTargetConProbAll_Depo[n1]) + type[n1].unitUtility_sur*(meetTargetConProbAll_Depo[n1]))*isContaminatedDepo;
					v_multicut[n1] -= b_ij * allocateBinVars[pools[j].depo_idx[i]][j];
				}
				meetTargetConProbAll_Depo.end();
			}
		}
		depo_prob_temp.end();
		meetTargetConProbAll_Pool.end();

		//actualUtility += poolUtility*poolNotContaminationProb;
		//for (n = 0; n < num_type; n++) 
		//	optCutRHS -= max_pool_size*(1 - poolNotContaminationProb)*type[n].unitUtility;
		//singlecutRHS -= max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
		//multicutRHS[j] = -max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
		//for (j1 = 0; j1<num_pool; j1++)
		//	v_multicut[j][j1] += max_pool_size*(1 - poolNotContaminationProb)*sumUnitUtility;
					//env.out() << v_multicut[j][0] << " " << multicutRHS[j] << endl;
	}

	for (n = 0; n < num_type; n++)
		v_multicut[n] -= expectedUtilityType[n];
	//env.out() << "actual utility" << actualUtility << ". expected utility" << expectedUtilityTypeSol[num_type] << endl;
	if (actualUtility_best < actualUtility) {
		actualUtility_best = actualUtility;
		for (i = 0; i < num_depo; i++) {
			for (j = 0; j < num_pool; j++) {
				for (n = 0; n < num_type; n++) {
					bestSolVol[i][j][n] = allocateVolVarsSol[i][j][n];
				}
			}
		}
	}
	if (actualUtility < expectedUtilityTypeSol[num_type] + RC_EPS) {
		violatedCutFound = IloTrue;
		//system("pause");
		cut_count++;
		if (cut_count % 10 == 0)
			env.out() << "Have added cuts ---- " << cut_count << endl;
	}
	for (j = 0; j < num_pool; j++)
		bigMArray[j].end();
	bigMArray.end();
	depoCleanProb.end();
	poolCleanProb.end();
	depoVol.end();
	productionTargets.end();
	meetTargetProbAll.end();


	//env.out() << actualUtility_best  << endl;
	//poolCleanProb.end();
	//productionTargets.end();
	//depoCleanProb.end();
	return violatedCutFound;

} // END OptCutAdded



void ObtainPoolInfo(const IloNumArray3 allocateVolVarsSol, vector<deposit> depo, vector<pool> &pools) {
	IloInt num_depo = depo.size();
	IloInt num_pool = pools.size();
	IloInt num_type = allocateVolVarsSol[0][0].getSize();
	//cout << "ObtainPoolInfo ===== " << num_depo << num_pool << num_type << endl;
	IloInt i, j, n;
	for (j = 0; j < num_pool; j++) {
		pools[j].is_used = false;
		pools[j].volume = 0;
		pools[j].depo_idx.clear();
		pools[j].clean_prob = 1;
		for (i = 0; i < num_depo; i++)
			for (n = 0; n < num_type; n++) {
				if (allocateVolVarsSol[i][j][n] > 1-RC_EPS) {
					//cout << "ObtainPoolInfo +++  " << i <<j<<n  << " " <<allocateVolVarsSol[i][j][n] <<endl;
					pools[j].is_used = true;
					pools[j].type_idx = n;
					pools[j].volume += allocateVolVarsSol[i][j][n];
					pools[j].depo_idx.push_back(i);
					pools[j].clean_prob *= depo[i].bacteriaFreeProb;
					//cout << "ObtainPoolInfo --- "<< j+1<<" "<<pools[j].is_used <<" "<< pools[j].type_idx<<" " <<pools[j].volume<<" "<<pools[j].clean_prob<<" "<<pools[j].depo_idx[pools[j].depo_idx.size()-1]<<endl;
					//	cout << "Pool " << j + 1 << " has " << allocateVolVarsSol[i][j][n] << " from deposit " << i + 1 << "producing type " << n + 1 << endl;
				}
			}
	}
}
/*
void UpdatePoolProb(vector<pool> &pools, IloNumArray depo_prob) {
	int num_pool = pools.size();
	for (int j = 0; j < num_pool; j++){
		if (pools[j].is_used){
			pools[j].clean_prob = 1;
			for (int i = 0; i < pools[j].depo_idx.size(); i++){
				pools[j].clean_prob *= (1 - depo_prob[pools[j].depo_idx[i]]);
			}
		}
	}
}

IloBoolArray ObtainTypeProduction(IloNumArray2 &TypePoolVol, IloNumArray2 &TypePoolProb, vector<pool> pools) {
	IloInt num_pool = pools.size();
	IloInt num_type = TypePoolVol.getSize();
	IloEnv env = TypePoolVol.getEnv();
	IloBoolArray nonEmptyType(env, num_type);
	IloInt n, j;
		for(n=0;n<num_type;n++ )
			nonEmptyType[n] = IloFalse;

	for (n = 0; n < num_type; n++) {
		TypePoolVol[n].clear();
		TypePoolProb[n].clear();
	}

	for (j = 0; j < num_pool; j++) {
		if(pools[j].is_used){
		TypePoolVol[pools[j].type_idx].add(pools[j].volume);
		TypePoolProb[pools[j].type_idx].add(pools[j].clean_prob);
		nonEmptyType[pools[j].type_idx] = IloTrue;
		}
	}
	return nonEmptyType;
}


IloNum ObtainMeetTargetProb(IloNumArray poolVol, IloNumArray poolProb, IloNum productionTarget) {
	if (productionTarget < RC_EPS) return 1;
	else if (poolVol.getSize() == 1) return (poolVol[0] < productionTarget) ? 0 : poolProb[0];
	else if (IloSum(poolVol) < productionTarget) return 0;
	else {
		IloNum currentPoolProb = poolProb[0];
		IloNum currentPoolVol = poolVol[0];
		poolProb.remove(0);
		poolVol.remove(0);
		return (1 - currentPoolProb)*ObtainMeetTargetProb(poolVol, poolProb, productionTarget) + currentPoolProb * ObtainMeetTargetProb(poolVol, poolProb, productionTarget - currentPoolVol);
	} 
}

IloNumArray ObtainMeetTargetProbAll(IloNumArray2 TypePoolVol, IloNumArray2 TypePoolProb, IloNumArray productionTarget, IloBoolArray producedType){
	IloInt num_type = TypePoolVol.getSize();
	IloEnv env = productionTarget.getEnv();
	IloNumArray meetTargetProb(env, num_type);
	IloInt n;
			for (int n = 0; n < num_type; n++){
			if (producedType[n])
				meetTargetProb[n] = ObtainMeetTargetProb(TypePoolVol[n] ,  TypePoolProb[n], productionTarget[n]);
			else if(productionTarget[n] < RC_EPS) meetTargetProb[n] = 1;
			else return 0;
		}
	return meetTargetProb;
}
*/

IloNumArray CalculateMeetTargetProb(vector<pool> pools, IloNumArray productionTargets, IloNumArray depoCleanProb, IloNumArray depoVol) {
	IloInt num_type = productionTargets.getSize();
	IloEnv env = depoCleanProb.getEnv();
	IloNum num_depo = depoCleanProb.getSize();
	IloNum num_pool = pools.size();
	//IloNumArray meetTargetProb(env, num_type);
	IloInt depo_idx = -1;
	IloNum largest_vol = 0;
	IloBoolArray depo_clean(env, num_depo);
	IloNum i;
	for (i = 0; i < num_depo; i++) {
		if (depoCleanProb[i]<1 - RC_EPS && depoCleanProb[i]> RC_EPS && depoVol[i] > largest_vol) { 
			depo_idx = i;
			largest_vol = depoVol[i];
		//	cout << " right before" << depo_idx << endl;
		}
		else depo_clean[i] = (depoCleanProb[i] > 1 - RC_EPS);
		//cout << i << " depo clean indicator " << depo_clean[i] << " " << depoCleanProb[i] << ". "<<endl;
	}
	if (depo_idx == -1) {
		IloNumArray result = WhetherMeetTarget( pools, productionTargets, depo_clean);
		depo_clean.end();
		return result;
	}
	else {
		IloNumArray depoCleanProb_copy1(env, num_depo);
		IloNumArray depoCleanProb_copy2(env, num_depo);
		for (i = 0; i < num_depo; i++) {
			depoCleanProb_copy1[i] = depoCleanProb[i];
			depoCleanProb_copy2[i] = depoCleanProb[i];
		}
		IloNum prob_temp = depoCleanProb[depo_idx];
		depoCleanProb_copy1[depo_idx] = 0;
		depoCleanProb_copy2[depo_idx] = 1;
		depo_clean.end();
		IloNumArray result2 = CalculateMeetTargetProb(pools, productionTargets, depoCleanProb_copy2, depoVol);
		IloNumArray result(env);
		if (IloSum(result2) < RC_EPS)
			return result2;
		else{
			result = SumofTwoArray(CalculateMeetTargetProb(pools, productionTargets, depoCleanProb_copy1, depoVol), 1- prob_temp, result2, prob_temp);
//			result1.end();
		}
		depoCleanProb_copy1.end();
		depoCleanProb_copy2.end();
		return result;
	}
}

IloNumArray WhetherMeetTarget(vector<pool> pools, IloNumArray productionTargets, IloBoolArray depo_clean) {
	IloInt num_type = productionTargets.getSize();
	IloEnv env = productionTargets.getEnv();
	IloNum num_depo = depo_clean.getSize();
	IloNum num_pool = pools.size();
	IloNum i, j, n;
	IloNumArray production(env, num_type);
	IloNumArray meetTarget(env, num_type);
	for (n = 0; n < num_type; n++) {
		production[n] = 0;
		meetTarget[n] = 1;
	}
	for (j = 0; j < num_pool; j++) {
		if (pools[j].is_used) {
			
			IloBool isclean=IloTrue;
			for (i = 0; i < pools[j].depo_idx.size(); i++){
				if (!depo_clean[pools[j].depo_idx[i]])
					isclean = IloFalse;
			}
			if (isclean) production[pools[j].type_idx] += pools[j].volume;
		}
	}
	for (n = 0; n < num_type; n++) {
		if (production[n] < productionTargets[n])
			meetTarget[n] = 0;
		//env.out() << " production[n] < productionTargets[n]" << production[n] << " " << productionTargets[n] << meetTarget[n] << endl;
	}
	production.end();
	//productionTargets.end();
	return meetTarget;
}

IloNumArray SumofTwoArray(IloNumArray array1, IloNum prob1, IloNumArray array2, IloNum prob2) {
	IloInt arraySize = array1.getSize();
	IloEnv env = array1.getEnv();
	IloNumArray result(env, arraySize);
	if (arraySize != array2.getSize()) {
		env.out() << arraySize << " " << prob1 << " " << array2.getSize() << " " << prob2 << endl;
		env.out() << "The sizes of two arrays don't match." << endl;
		return result;
	}
	for (IloInt i = 0; i < arraySize; i++) {
		result[i] = array1[i] * prob1 + array2[i] * prob2;
	}
	array1.end();
	array2.end();
	return result;
}

IloNum evaluateActualUtility(IloNumArray depoCleanProb, vector<pool> pools,vector<milktype> type) {
	IloInt num_type = type.size();
	IloEnv env = depoCleanProb.getEnv();
	IloNum num_pool = pools.size();
	IloNum num_depo = depoCleanProb.getSize();
	IloInt depo_idx = -1;
	IloNum i, j, n;
	for (i = 0; i < num_depo; i++) {
		if (depoCleanProb[i]<1 - RC_EPS && depoCleanProb[i]> RC_EPS)
			depo_idx = i;
			//	cout << " right before" << depo_idx << endl;
		//cout << i << " depo clean indicator " << depo_clean[i] << " " << depoCleanProb[i] << ". "<<endl;
	}
	if (depo_idx == -1) {
		IloNum actualUtility = 0;
		IloNumArray production(env, num_type);
		for (n = 0; n < num_type; n++) {
			production[n] = 0;
		}
		for (j = 0; j < num_pool; j++) {
			if(pools[j].is_used){
				IloBool poolIsContaminated = IloFalse;
				for (int i = 0; i < pools[j].depo_idx.size(); i++) {
					if (depoCleanProb[pools[j].depo_idx[i]] < RC_EPS) 
						poolIsContaminated = IloTrue;
				}
				if(!poolIsContaminated)
					production[pools[j].type_idx] += pools[j].volume;
		}
		}
		for (n = 0; n < num_type; n++) {
			if (production[n] > type[n].productionTarget)
				actualUtility += (production[n] - type[n].productionTarget)*type[n].unitUtility_sur + type[n].productionTarget*type[n].unitUtility_def;
			else actualUtility += production[n] * type[n].unitUtility_def;
		}
		//env.out() << depoCleanProb <<" actualUtility " << actualUtility << endl;
		production.end();
		return actualUtility;
	}
	else {
		IloNumArray depoCleanProb_copy1(env, num_depo);
		IloNumArray depoCleanProb_copy2(env, num_depo);
		for (i = 0; i < num_depo; i++) {
			depoCleanProb_copy1[i] = depoCleanProb[i];
			depoCleanProb_copy2[i] = depoCleanProb[i];
		}
		//env.out() << depo_idx << " ";
		IloNum prob_temp = depoCleanProb[depo_idx];
		depoCleanProb_copy1[depo_idx] = 0;
		depoCleanProb_copy2[depo_idx] = 1;
		IloNum result2 = evaluateActualUtility(depoCleanProb_copy2, pools, type);
		IloNum result;
		if (result2 < RC_EPS)
			return 0;
		//env.out() << (1 - prob_temp) <<" "<< depoCleanProb_copy1 << " evaluateActualUtility " << evaluateActualUtility(depoCleanProb_copy1, pools, type) << endl;
		else 
			result = (1 - prob_temp)*evaluateActualUtility(depoCleanProb_copy1, pools, type) + prob_temp* result2;
		depoCleanProb_copy1.end();
		depoCleanProb_copy2.end();
		return result;
	}
}

double evaluateUtilityBySol(vector<milktype> type, vector<double> depoCleanProb, double *** Sol) {
	int num_type = type.size();
	int num_depo = depoCleanProb.size();
	int num_pool = num_depo;
	int i, j, n;
	int depo_idx =-1;

	for (i = 0; i < num_depo; i++) {
		if (depoCleanProb[i]<1 - RC_EPS && depoCleanProb[i]> RC_EPS)
			depo_idx = i;
		//	cout << " right before" << depo_idx << endl;
	//cout << i << " depo clean indicator " << depo_clean[i] << " " << depoCleanProb[i] << ". "<<endl;
	}

	if (depo_idx == -1) {
		int * productionByType = new int[num_type];
		for (n = 0; n < num_type; n++)
			productionByType[n] = 0;

		for (j = 0; j < num_pool; j++) {
			int poolProduction = 0;
			int typeDetermined = -1;
			for (n = 0; n < num_pool; n++) {
				//			bool poolContaminated = false;
				for (i = 0; i < num_depo; i++) {
					if (Sol[i][j][n] > RC_EPS) {
						typeDetermined = n;
						if (depoCleanProb[i] < RC_EPS) {
							//							poolContaminated = true;

							poolProduction = 0;
							break;
						}
						else poolProduction += Sol[i][j][n];
					}
				}
				if (typeDetermined > -1) break;
			}
			if (typeDetermined > -1)
				productionByType[typeDetermined] += poolProduction;
		}
		double utility = 0;
		for (n = 0; n < num_type; n++) {
			if (productionByType[n] > type[n].productionTarget)
				utility += (productionByType[n] - type[n].productionTarget)*type[n].unitUtility_sur + type[n].productionTarget*type[n].unitUtility_def;
			else
				utility += productionByType[n] * type[n].unitUtility_def;
		}
		delete[] productionByType;
		return utility;
	}
	else {
		vector<double> depoCleanProb_copy1(depoCleanProb);
		vector<double> depoCleanProb_copy2(depoCleanProb);
		//for (i = 0; i < num_depo; i++) {
		//	depoCleanProb_copy1.push_back(depoCleanProb[i]);
		//	depoCleanProb_copy2.push_back(depoCleanProb[i]);
		//}
		double prob_temp = depoCleanProb[depo_idx];
		depoCleanProb_copy1[depo_idx] = 1;
		depoCleanProb_copy2[depo_idx] = 0;
		double utility = prob_temp * evaluateUtilityBySol(type, depoCleanProb_copy1, Sol) + (1 - prob_temp)*evaluateUtilityBySol(type, depoCleanProb_copy2, Sol);
		depoCleanProb_copy1.clear();
		depoCleanProb_copy2.clear();
		return utility;
	}
}
