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
const string MFSEDTABLE = "../lookupTables/MFSED_nB_3_dQ_0.001_dimF_512.txt";
const string FBORDERS = "../lookupTables/Borders_nB_3_dimF_512.txt";
const size_t NROWS = 8;
const int NEVALKEY = 8; 
const int NFEAT = 512;
const int64_t THRESHOLD = 2000;
const int SUBJECTS = 500;
const size_t SAMPLES = 3;
const int SEED = 45676;

int main(int argc, char *argv[]) {
    TimeVar runTemp, runProbe, runVerif;

    string sBits = argv[1];

    string RUNTIMEFILE = "../experiments/runtimesMFSED/MFSEDClearComp/expMFSED-ClearComp"+sBits+".csv";
    string TEMPRUNTIMEFILE = "../experiments/runtimesMFSED/MFSEDClearComp/expMFSED-ClearComp-Template"+sBits+"-NEVALKEY-"+to_string(NEVALKEY)+".csv";
    string PROBRUNTIMEFILE = "../experiments/runtimesMFSED/MFSEDClearComp/expMFSED-ClearComp-Probe"+sBits+".csv";

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



    // Generate the rotation evaluation keys
    auto shiftIndexes = genVectOfInt(1, NEVALKEY);
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, shiftIndexes);


    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    cout << "ringDim = " << ringDim << endl;


    // Loading borders and the MFIP lookup table
    auto borders = readLineDouble(FBORDERS);
    auto tabMFSED = readTableInt64(MFSEDTABLE, NROWS);

    size_t nPackedRowsPerCipher = (size_t) NEVALKEY/NROWS;
    size_t numCiphers = (size_t) (NFEAT * NROWS)/NEVALKEY;

    cout << "numCiphers = " << numCiphers << "; nPackedRowPerCipher = " << nPackedRowsPerCipher << endl;




    vector<Ciphertext<DCRTPoly>> subj1RefTempCt, subj1ProbeCt, subj2ProbeCt;   



    for (int subject1 = 1; subject1 < SUBJECTS+1; subject1++) {

        
         
        auto permutations = genPermutations(SEED, numCiphers, NEVALKEY);
        string subj1RefTempF = SAMPLESDIR+to_string(subject1)+"_1.txt";
        string subj1ProbeF = SAMPLESDIR+to_string(subject1)+"_2.txt";

        cout << "subj1RefTempF = " << subj1RefTempF << endl;
        cout << "subj1ProbeF = " << subj1ProbeF << endl;

        

        /////////////////////
        // Enrollment Phase
        /////////////////////      
        
        // TIC(runTemp);                
        auto subj1RefTempS = readLineDouble(subj1RefTempF);
        auto subj1RefTempQ = quantizeFeatures(borders, subj1RefTempS);
        auto subj1RefTempPacked = packRefTemplate2(subj1RefTempQ, tabMFSED, nPackedRowsPerCipher, numCiphers);        
        subj1RefTempCt = genEncryptedReferenceTemplate(cryptoContext, keyPair, subj1RefTempPacked, permutations);
        // writeProcTime(TEMPRUNTIMEFILE, TOC_MS(runTemp));




        // For isolating the first plaintext slot
        vector<int64_t> vectR = genRandVectOfInt64(1, (int64_t) plaintextModulus/4-1,ringDim);
        vectR.at(0) = 0;
        Plaintext blindRPt = cryptoContext->MakePackedPlaintext(vectR);
        auto blindRCt = cryptoContext->Encrypt(keyPair.publicKey, blindRPt);



        /////////////////////
        // Verification Phase
        /////////////////////  

        auto subj1ProbeS = readLineDouble(subj1ProbeF);
        auto subj1ProbeQ = quantizeFeatures(borders, subj1ProbeS);



        // TIC(runProbe);
        auto permProbInd = genPackedProbeIndexes(subj1ProbeQ, permutations, nPackedRowsPerCipher, NROWS);
        writeProcTime(PROBRUNTIMEFILE, TOC_MS(runProbe));
        


        // Clear-text comparison with the threshold (Dissimilarity Score Test: S <= threshold ?)
        
        TIC(runVerif);
        auto finalScoreCt = computeFinalScoreMFBRClearComp(cryptoContext, subj1RefTempCt, permProbInd, blindRCt, nPackedRowsPerCipher, NROWS);
        Plaintext finalScorePT;
        cryptoContext->Decrypt(keyPair.secretKey, finalScoreCt, &finalScorePT);
        auto fs = finalScorePT->GetPackedValue().at(0);
        string resultStr = (THRESHOLD<fs)? "No Match" : "Match"; 
        auto t = TOC_MS(runVerif); 
        writeProcTime(RUNTIMEFILE, t);
        cout << resultStr << " => Final Dissimilarity Score = " << fs << endl;  
        cout << "Verification Time : " << t << "ms" << endl; 



    }


    return 0;
}