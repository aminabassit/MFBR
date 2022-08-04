#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>


#include "palisade.h"
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "pubkeylp-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"

#include "include/mfbr.hpp"

using namespace std;
using namespace lbcrypto;

const string SAMPLESDIR = "../data/VGGFace2/ArcFace-R100/ArcFace-R100_";
const string FBORDERS = "../lookupTables/Borders_nB_3_dimF_512.txt";
const int NFEAT = 512;
const int SEED = 45676;



int main(int argc, char *argv[]) {

    string sBits = argv[1];
    string btpProbeFolder = "../experiments/templates/probe/IPandSED"+sBits+"/";
    string btpReferenceFolder = "../experiments/templates/reference/IPandSED"+sBits+"/";
    string btpProbeFile("Probe-");
    string btpReferenceFile("Reference-");
    string subj1RefTempF = SAMPLESDIR+"1_1.txt";    
    string subj1ProbeF = SAMPLESDIR+"1_2.txt";

    cout << "HEStd_"+sBits+"_classic" << endl;

    
    int64_t nBits(16);
    uint64_t m(32768);
    auto plaintextMod = FirstPrime<NativeInteger>(nBits, m);
    usint plaintextModulus = (usint) plaintextMod.ConvertToInt();
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
        numMults = 2;
        numAdds = 0;
        numKeyswitches = 0;
        dcrtBits = 36; 
    }
    if (sBits == "192")
    {
        securityLevel = HEStd_192_classic;
        numMults = 2;
        numAdds = 0;
        numKeyswitches = 0;
        dcrtBits = 37; 
    }
    if (sBits == "256")
    {
        securityLevel = HEStd_256_classic;
        numMults = 2;
        numAdds = 0;
        numKeyswitches = 0;
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


    auto shiftIndexes = genVectOfInt(1, NFEAT);
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, shiftIndexes);


    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    cout << "ringDim = " << ringDim << endl;

    // Loading 
    auto borders = readLineDouble(FBORDERS);

    



    auto subj1RefTempS = readLineDouble(subj1RefTempF);
    auto subj1RefTempQ = quantizeFeatures(borders, subj1RefTempS);        
    auto subj1RefTempCt = encryptVectorInt64(cryptoContext, keyPair, subj1RefTempQ);
    cout << "Template quantized and encrypted..." << endl;

    if (!Serial::SerializeToFile(btpReferenceFolder + btpReferenceFile + ".txt", subj1RefTempCt, SerType::BINARY)) {
        std::cerr
            << "Error writing serialization of a baseline reference to "<< btpReferenceFile <<".txt"
            << std::endl;
        return 1;
    } 

    cout << "Baseline reference saved ......." << endl; 


    // Probe generation
    
    auto subj1ProbeS = readLineDouble(subj1ProbeF);
    auto subj1ProbeQ = quantizeFeatures(borders, subj1ProbeS);
    auto subj1ProbeCT = encryptVectorInt64(cryptoContext, keyPair, subj1ProbeQ);
    cout << "Probe quantized and encrypted..." << endl;


    if (!Serial::SerializeToFile(btpProbeFolder + btpProbeFile + ".txt", subj1ProbeCT, SerType::BINARY)) {
        std::cerr
            << "Error writing serialization of a baseline probe to "<< btpProbeFile <<".txt"
            << std::endl;
        return 1;
    } 

    cout << "Baseline probe saved ......." << endl; 



    return 0;
}









