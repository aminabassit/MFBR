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
const string FBORDERS = "../lookupTables/Borders_nB_3_dimF_512.txt";
const int NFEAT = 512;
const int64_t THRESHOLD = 1;
const int SUBJECTS = 500;
const size_t SAMPLES = 3;



int main(int argc, char *argv[]) {
    TimeVar runTemp, runProbe, runVerif;

    string sBits = argv[1];
    string RUNTIMEFILE = "../experiments/runtimesIPBoddeti/IPBoddetiClearComp/expIPBoddeti-ClearComp"+sBits+".csv";
    string TEMPRUNTIMEFILE = "../experiments/runtimesBoddeti/IPBoddetiClearComp/expIPBoddeti-ClearComp-Template"+sBits+".csv";
    string PROBRUNTIMEFILE = "../experiments/runtimesBoddeti/IPBoddetiClearComp/expIPBoddeti-ClearComp-Probe"+sBits+".csv";

    cout << "HEStd_"+sBits+"_classic" << endl;
    
    int precision = 400; // 0.0025;

    int64_t nBits(20); 
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
        numMults = 1;
        numAdds = 0;
        numKeyswitches = 0;
        dcrtBits = 36; 
    }
    if (sBits == "192")
    {
        securityLevel = HEStd_192_classic;
        numMults = 1;
        numAdds = 0;
        numKeyswitches = 0;
        dcrtBits = 37; 
    }
    if (sBits == "256")
    {
        securityLevel = HEStd_256_classic;
        numMults = 1;
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

    // Generate the rotation evaluation keys
    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    cout << "ringDim = " << ringDim << endl;

    auto shiftIndexes = genVectOfPowOf2(log2(ringDim)+1);
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, shiftIndexes);

    size_t row_size = (size_t) ringDim/2;

    vector<Ciphertext<DCRTPoly>> subj1RefTempCt, subj1ProbeCt;  



    for (int subject1 = 1; subject1 < SUBJECTS+1; subject1++) {

        
        string subj1RefTempF = SAMPLESDIR+to_string(subject1)+"_1.txt";
        string subj1ProbeF = SAMPLESDIR+to_string(subject1)+"_2.txt";

        cout << "subj1RefTempF = " << subj1RefTempF << endl;
        cout << "subj1ProbeF = " << subj1ProbeF << endl;

        

        /////////////////////
        // Enrollment Phase
        /////////////////////
 

        // TIC(runTemp);
        auto subj1RefTempS = readLineDouble(subj1RefTempF);
        auto subj1RefTempQ = quantizeFeaturesBoddeti(precision, subj1RefTempS);      
        auto subj1RefTempCt = encryptVectorInt64(cryptoContext, keyPair, subj1RefTempQ);
        // writeProcTime(TEMPRUNTIMEFILE, TOC_MS(runTemp));



        /////////////////////
        // Verification Phase
        /////////////////////
        
        
        auto subj1ProbeS = readLineDouble(subj1ProbeF);
        auto subj1ProbeQ = quantizeFeaturesBoddeti(precision, subj1ProbeS);


        
        // TIC(runProbe);
        auto subj1ProbeCT = encryptVectorInt64(cryptoContext, keyPair, subj1ProbeQ);
        // writeProcTime(PROBRUNTIMEFILE, TOC_MS(runProbe));
        

        // Clear-text comparison with the threshold (Similarity Score Test: threshold <= S ?)

        TIC(runVerif);
        auto finalScoreCt = computeFinalScoreIPBoddeti(cryptoContext, subj1RefTempCt, subj1ProbeCT, row_size);
        Plaintext finalScorePT;
        cryptoContext->Decrypt(keyPair.secretKey, finalScoreCt, &finalScorePT);
        // printVectorInt64(finalScorePT->GetPackedValue());
        auto fs = finalScorePT->GetPackedValue().at(0);
        string resultStr = (fs<THRESHOLD)? "No Match" : "Match";       
        auto t = TOC_MS(runVerif); 
        writeProcTime(RUNTIMEFILE, t);
        cout << resultStr << " => Final Similarity Score = " << fs << endl; 
        cout << "Verification Time : " << t << "ms" << endl;   


    }





    return 0;
}