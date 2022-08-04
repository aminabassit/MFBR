#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>


#include "palisade.h"

#include "include/mfbr.hpp"

using namespace std;
using namespace lbcrypto;


const string FBORDERS = "../lookupTables/Borders_nB_3_dimF_512.txt";
const string SAMPLESDIR = "../data/VGGFace2/ArcFace-R100/ArcFace-R100_";
const int NFEAT = 512;
const int64_t THRESHOLD = 5498; 
const int SUBJECTS = 500;
const size_t SAMPLES = 3;



int main(int argc, char *argv[]) {
    TimeVar runTemp, runProbe, runVerif;

    string sBits = argv[1];
    string RUNTIMEFILE = "../experiments/runtimesBaseline/IPBaselineClearComp/expIP-ClearComp"+sBits+".csv";
    string TEMPRUNTIMEFILE = "../experiments/runtimesBaseline/IPBaselineClearComp/expIP-ClearComp-Template"+sBits+".csv";
    string PROBRUNTIMEFILE = "../experiments/runtimesBaseline/IPBaselineClearComp/expIP-ClearComp-Probe"+sBits+".csv";

    cout << "HEStd_"+sBits+"_classic" << endl;

    
    int64_t nBits(16); 
    uint64_t m(32768);
    auto plaintextMod = FirstPrime<NativeInteger>(nBits, m);
    usint plaintextModulus = (usint) plaintextMod.ConvertToInt();
    std::cout << "plaintextModulus = " << plaintextModulus << std::endl;

    EncodingParams encodingParams(std::make_shared<EncodingParamsImpl>(plaintextModulus));
    double sigma = 3.2;
    usint numMults(0), numAdds(0), numKeyswitches(0);
    int maxDepth(2);
    uint32_t relinWindow(0);
    size_t dcrtBits(30);
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

    
    // Generate the rotation evaluation keys
    auto shiftIndexes = genVectOfInt(1, NFEAT);
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, shiftIndexes);


    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    cout << "ringDim = " << ringDim << endl;


    // Loading borders
    auto borders = readLineDouble(FBORDERS);    
    




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
        auto subj1RefTempQ = quantizeFeatures(borders, subj1RefTempS);              
        auto subj1RefTempCt = encryptVectorInt64(cryptoContext, keyPair, subj1RefTempQ);
        // writeProcTime(TEMPRUNTIMEFILE, TOC_MS(runTemp));
        


        
        // For isolating the first plaintext slot 
        vector<int64_t> vectOneAndZeros(ringDim, 0); 
        vectOneAndZeros.at(0) = 1;
        auto oneAndZerosPT = cryptoContext->MakePackedPlaintext(vectOneAndZeros);
        auto oneAndZerosCT = cryptoContext->Encrypt(keyPair.publicKey, oneAndZerosPT);
        


        /////////////////////
        // Verification Phase
        /////////////////////

        
        auto subj1ProbeS = readLineDouble(subj1ProbeF);
        auto subj1ProbeQ = quantizeFeatures(borders, subj1ProbeS); 

        
        
        // TIC(runProbe);
        auto subj1ProbeCT = encryptVectorInt64(cryptoContext, keyPair, subj1ProbeQ);
        // writeProcTime(PROBRUNTIMEFILE, TOC_MS(runProbe));       


        // Clear-text comparison with the threshold (Similarity Score Test: threshold <= S ?)

        TIC(runVerif);
        auto finalScoreCt = computeFinalScoreIPBaseline(cryptoContext, subj1RefTempCt, subj1ProbeCT, oneAndZerosCT, NFEAT);
        Plaintext finalScorePT;
        cryptoContext->Decrypt(keyPair.secretKey, finalScoreCt, &finalScorePT);
        auto fs = finalScorePT->GetPackedValue().at(0);
        string resultStr = (fs<THRESHOLD)? "No Match" : "Match";       
        auto t = TOC_MS(runVerif); 
        writeProcTime(RUNTIMEFILE, t);
        cout << resultStr << " => Final Similarity Score = " << fs << endl; 
        cout << "Verification Time : " << t << "ms" << endl;      
        

    }





    return 0;
}