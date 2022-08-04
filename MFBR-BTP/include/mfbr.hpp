#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <random>
#include <iomanip>
#include <list>
#include <numeric>






#include "cryptocontextgen.h"
#include "palisade.h"

using namespace std;

void writeProcTime(string filename, double processingTime);
const vector<string> explode(const string& s, const char& c);
vector<int> genVectOfInt(int begin, size_t len);
vector<int> genVectOfPowOf2(size_t len);
vector<double> readLineDouble(string filename);
map< int, vector<int64_t>> readTableInt64(string filename, size_t nRows);
vector<int64_t> permuteSample(vector<int64_t> sample, vector<int> permutation);
map<int, vector<int>> genPermutations(int seed, size_t nPerm, size_t lenPerm);
vector<int64_t> vectIntPermutationInt64Neg(int seed, int64_t begin, size_t len);
vector<int64_t> genVectOfInt64Neg(int64_t begin, size_t len);
vector<int> vectIntPermutation(int seed, int begin, size_t len);
vector<vector<int>> getIndexRotationGroups(vector<int> permProbInd);




vector<int64_t> quantizeFeatures(vector<double> borders, vector<double> unQuantizedFeatures);
vector<int64_t> quantizeFeaturesBoddeti(int precision, vector<double> unQuantizedFeatures);


Ciphertext<DCRTPoly> getRotationGroupAtIndex(CryptoContext<DCRTPoly> cryptoContext, vector<Ciphertext<DCRTPoly>> encRefTempCT, vector<int> rotaionIndexesInd, int ind);
Ciphertext<DCRTPoly> encryptVectorInt64(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, vector<int64_t> vect);


vector<int64_t> genRandVectOfInt64(int64_t begin, int64_t end, size_t len);
vector<vector<int64_t>> packRefTemplate2(vector<int64_t> sampleQ, map<int, vector<int64_t>> mfbrTab, size_t nPackedRowsPerCipher, size_t numCiphers);
vector<int> genPackedProbeIndexes(vector<int64_t> probSampleQ, map<int, vector<int>> permutations, size_t nPackedRowsPerCipher, size_t nRows);
vector<Ciphertext<DCRTPoly>> genEncryptedReferenceTemplate(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, vector<vector<int64_t>> reftemp, map<int, vector<int>> permutations);



Ciphertext<DCRTPoly> computeFinalScoreIPBaseline(CryptoContext<DCRTPoly> cryptoContext, Ciphertext<DCRTPoly> encRefTempCT, Ciphertext<DCRTPoly> encProbeCT, Ciphertext<DCRTPoly> oneAndZerosCT, size_t nFeat);
Ciphertext<DCRTPoly> computeFinalScoreSEDBaseline(CryptoContext<DCRTPoly> cryptoContext, Ciphertext<DCRTPoly> encRefTempCT, Ciphertext<DCRTPoly> encProbeCT, Ciphertext<DCRTPoly> oneAndZerosCT, size_t nFeat);

Ciphertext<DCRTPoly> computeFinalScoreMFBRClearComp(CryptoContext<DCRTPoly> cryptoContext, vector<Ciphertext<DCRTPoly>> encRefTempCT, vector<int> permProbInd, Ciphertext<DCRTPoly> blindRCT, size_t nPackedRowsPerCipher, int nRows);
Ciphertext<DCRTPoly> computeFinalScoreForBlindedCompGroups(CryptoContext<DCRTPoly> cryptoContext, vector<Ciphertext<DCRTPoly>> encRefTempCT, vector<int> permProbInd, Plaintext onePlain);

Ciphertext<DCRTPoly> computeFinalScoreIPBoddeti(CryptoContext<DCRTPoly> cryptoContext, Ciphertext<DCRTPoly> encRefTempCT, Ciphertext<DCRTPoly> encProbeCT, size_t row_size);






bool verifyComparisonInt64(vector<int64_t> blindedCompVect);
bool verifyComparison(Plaintext blindedCompVectPlain, size_t len);




