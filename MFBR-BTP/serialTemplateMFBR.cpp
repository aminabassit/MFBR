#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "cryptocontextgen.h"
#include "palisade.h"
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "pubkeylp-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"

#include "include/mfbr.hpp"

using namespace std;
using namespace lbcrypto;


const string SAMPLESDIR = "../data/VGGFace2/ArcFace-R100/ArcFace-R100_";
const string MFIPTABLE = "../lookupTables/MFIP_nB_3_dQ_0.001_dimF_512.txt";
const string FBORDERS = "../lookupTables/Borders_nB_3_dimF_512.txt";
const size_t NROWS = 8;
const int NEVALKEY = 8; 
const int NFEAT = 512;
const int64_t THRESHOLD = 131;
const int SUBJECTS = 500;
const size_t SAMPLES = 3;
const int SEED = 45676;


int main(int argc, char *argv[]) {

    
    string sBits = argv[1];
    string btpFolder = "../experiments/templates/reference/MFBR"+sBits+"/";
    string btpFile("Reference-");


    cout << "HEStd_"+sBits+"_classic" << endl;

    int64_t nBits(16);
    uint64_t m(32768);
    auto plaintextMod = FirstPrime<NativeInteger>(nBits, m);
    int plaintextModulus = (int) plaintextMod.ConvertToInt();
    std::cout << "plaintextModulus = " << plaintextModulus << std::endl;
    EncodingParams encodingParams(std::make_shared<EncodingParamsImpl>(plaintextModulus));
    double sigma = 3.2;
    usint numMults, numAdds, numKeyswitches;
    int maxDepth(2);
    uint32_t relinWindow(0);
    size_t dcrtBits;
    uint32_t n(0);
    SecurityLevel securityLevel;




    if (sBits == "128")
    {
        securityLevel = HEStd_128_classic;
        numMults = 0;
        numAdds = 0;
        numKeyswitches = 1;
        dcrtBits = 36; 
    }
    if (sBits == "192")
    {
        securityLevel = HEStd_192_classic;
        numMults = 0;
        numAdds = 0;
        numKeyswitches = 1;
        dcrtBits = 37; 
    }
    if (sBits == "256")
    {
        securityLevel = HEStd_256_classic;
        numMults = 0;
        numAdds = 0;
        numKeyswitches = 1;
        dcrtBits = 38; 
    } 


    CryptoContext<DCRTPoly> cryptoContext =
    CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
            encodingParams, securityLevel, sigma, numAdds, numMults, numKeyswitches, OPTIMIZED, maxDepth, relinWindow, dcrtBits, n);



    


    


    cryptoContext->Enable(ENCRYPTION);
    cryptoContext->Enable(SHE);


    LPKeyPair<DCRTPoly> keyPair;
    keyPair = cryptoContext->KeyGen();

    cryptoContext->EvalMultKeysGen(keyPair.secretKey);
    cryptoContext->EvalSumKeyGen(keyPair.secretKey);    

    // // Generate the rotation evaluation keys


    auto shiftIndexes = genVectOfInt(1, NEVALKEY);
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, shiftIndexes);


    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    cout << "ringDim = " << ringDim << endl;


    // Loading borders and the MFIP lookup table 
    auto borders = readLineDouble(FBORDERS);
    auto tabMFIP = readTableInt64(MFIPTABLE, NROWS);

    size_t nPackedRowsPerCipher = (size_t) NEVALKEY/NROWS;
    size_t numCiphers = (size_t) (NFEAT * NROWS)/NEVALKEY;

    cout << "numCiphers = " << numCiphers << "; nPackedRowPerCipher = " << nPackedRowsPerCipher << endl;
    
    cout << "Table and Borders were loaded..." << endl;

    auto permutations = genPermutations(SEED, numCiphers, NEVALKEY);
    string subj1RefTempF = SAMPLESDIR+"1_1.txt";

    cout << "subj1RefTempF = " << subj1RefTempF << endl;

    

    // Template encryption         
    
    auto subj1RefTempS = readLineDouble(subj1RefTempF);
    auto subj1RefTempQ = quantizeFeatures(borders, subj1RefTempS);
    auto subj1RefTempPacked = packRefTemplate2(subj1RefTempQ, tabMFIP, nPackedRowsPerCipher, numCiphers);
    cout << "Template packed..." << endl;
    auto subj1RefTempCt = genEncryptedReferenceTemplate(cryptoContext, keyPair, subj1RefTempPacked, permutations);





    for (size_t i = 0; i < subj1RefTempCt.size(); i++)
    {
        if (!Serial::SerializeToFile(btpFolder + btpFile + to_string(i) + ".txt", subj1RefTempCt.at(i), SerType::BINARY)) {
        std::cerr
            << "Error writing serialization of a MFBR reference template to "<< btpFile <<".txt"
            << std::endl;
        return 1;
        }   
    }
        
    cout << "MFBR reference template saved ......." << endl; 

    return 0;
}