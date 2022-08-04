#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>


#include "palisade.h"

#include "include/mfbr.hpp"

using namespace std;
using namespace lbcrypto;

const string SAMPLESDIR = "../data/VGGFace2/ArcFace-R100/ArcFace-R100_";
const string MFIPTABLE = "../lookupTables/MFIP_nB_3_dQ_0.001_dimF_512.txt";
const string FBORDERS = "../lookupTables/Borders_nB_3_dimF_512.txt";
const size_t NROWS = 8;
const int NEVALKEY = 8; 
const int NFEAT = 512;
const int SUBJECTS = 500;
const size_t SAMPLES = 3;
const int SEED = 45676;



int main(int argc, char *argv[]) {

    
    string btpFolder = "../experiments/templates/probe/MFBR/";
    string btpFile("Probe-");

    // Loading 
    auto borders = readLineDouble(FBORDERS);
    auto tabMFIP = readTableInt64(MFIPTABLE, NROWS);

    size_t nPackedRowsPerCipher = (size_t) NEVALKEY/NROWS;
    size_t numCiphers = (size_t) (NFEAT * NROWS)/NEVALKEY;


    auto permutations = genPermutations(SEED, numCiphers, NEVALKEY);
    string subj1ProbeF = SAMPLESDIR+"1_2.txt";

    auto subj1ProbeS = readLineDouble(subj1ProbeF);
    auto subj1ProbeQ = quantizeFeatures(borders, subj1ProbeS);
    cout << "Probe quantized..." << endl;
    auto permProbInd = genPackedProbeIndexes(subj1ProbeQ, permutations, nPackedRowsPerCipher, NROWS);


    for (size_t i = 0; i < permProbInd.size(); i++)
    {
        if (!Serial::SerializeToFile(btpFolder + btpFile + to_string(i) + ".txt", permProbInd[i], SerType::BINARY)) {
        std::cerr
            << "Error writing serialization of a MFBR probe to "<< btpFile <<".txt"
            << std::endl;
        return 1;
        }   
    }




    cout << "MFBR probe saved ......." << endl; 



    return 0;
}









