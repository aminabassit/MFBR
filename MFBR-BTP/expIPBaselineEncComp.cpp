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
const size_t COMPLEN = 4146; // number of scores between the threshold to Smax
const int64_t THRESHOLD = 5498;
const int SUBJECTS = 500;
const size_t SAMPLES = 3;
const int SEED = 45676;


int main(int argc, char *argv[]) {
    TimeVar runTemp, runProbe, runVerif;

    string sBits = argv[1];
    string RUNTIMEFILE = "../experiments/runtimesBaseline/IPBaselineEncComp/expIP-EncComp"+sBits+".csv";
    string TEMPRUNTIMEFILE = "../experiments/runtimesBaseline/IPBaselineEncComp/expIP-EncComp-Template"+sBits+".csv";
    string PROBRUNTIMEFILE = "../experiments/runtimesBaseline/IPBaselineEncComp/expIP-EncComp-Probe"+sBits+".csv";

    cout << "HEStd_"+sBits+"_classic" << endl;

    

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
        numMults = 2;
        numAdds = 0;
        numKeyswitches = 0;
        dcrtBits = 37; 
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

    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    cout << "ringDim = " << ringDim << endl;
    size_t row_size = (size_t) ringDim/2;   

    // Generate the rotation evaluation keys
    auto pow2Vect = genVectOfPowOf2(log2(row_size)+1);
    auto shiftIndexes = genVectOfInt(1, NFEAT);
    shiftIndexes.insert(shiftIndexes.end(), pow2Vect.begin(), pow2Vect.end());
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, shiftIndexes);


    


    // Loading borders
    auto borders = readLineDouble(FBORDERS);


    

    // For the encrypted comparison with the threshold: generate the comparison vector permuted and the blinding vector
    auto vectR = genRandVectOfInt64(1, (int64_t) plaintextModulus/4-1, row_size);
    Plaintext blindRPt = cryptoContext->MakePackedPlaintext(vectR);
    auto blindRCT = cryptoContext->Encrypt(keyPair.publicKey, blindRPt);
    vector<Ciphertext<DCRTPoly>> compVectTCs;    
    for (size_t i = 0; i < COMPLEN; i = i+row_size)
    {
        auto compVect = vectIntPermutationInt64Neg(SEED, THRESHOLD+(int64_t)i, row_size);
        Plaintext compVectPT = cryptoContext->MakePackedPlaintext(compVect);
        auto compVectCT = cryptoContext->Encrypt(keyPair.publicKey, compVectPT);
        compVectTCs.push_back(compVectCT);
    }
    cout << "len compVectTCs = " << compVectTCs.size() << endl;




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



        /**
         * Encrypted comparison with the threshold following the proceeder described in [BHP+]
         * Note that this proceeder tests a similarity score 
         * Test: threshold <= S ?)
         * TODO: Encrypted comparison with the threshold for testing a dissimilarity score (for SED)
         * [BHP+] https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9585508
        */
        


        TIC(runVerif);
        auto finalScoreCt0 = computeFinalScoreIPBaseline(cryptoContext, subj1RefTempCt, subj1ProbeCT, oneAndZerosCT, NFEAT);  

        // replicate the final score 'S' over all plaintext slots to get (S, S, ..., S)
        Ciphertext<DCRTPoly> finalScoreCT1 = finalScoreCt0;       
        for (size_t i = 0; i < log2(row_size); i++)
        {
            finalScoreCt0 = cryptoContext->EvalAtIndex(finalScoreCT1, pow(2,i)); 
            cryptoContext->EvalAddInPlace(finalScoreCT1, finalScoreCt0);
        }

        // compute the blinded comparison vector (.., (S-theta_{i})*r_{i}, ...)_{i} 
        // where theta_{i} \in [T, Smax] and r_{i} is rand used for the blinding 
        size_t lenCom = compVectTCs.size();
        vector<Ciphertext<DCRTPoly>> finalScoreBlindedCT(lenCom);
        #pragma omp parallel for schedule(dynamic) 
        for (size_t i = 0; i < lenCom; i++)               
        {
            auto finalScoreCT = cryptoContext->EvalAdd(finalScoreCT1, compVectTCs.at(i));
            finalScoreBlindedCT.at(i) = cryptoContext->EvalMult(finalScoreCT, blindRCT);  
        }

        bool res = false;

        for (size_t i = 0; i < lenCom; i++)
        {
            Plaintext finalScorePT;
            cryptoContext->Decrypt(keyPair.secretKey, finalScoreBlindedCT.at(i), &finalScorePT);
            finalScorePT->SetLength(row_size);
            auto fsc = finalScorePT->GetPackedValue();
            res = verifyComparisonInt64(fsc);
            if (res == true)
            {
                break;
            }
            
        }      
        
        string resultStr = (res)? "Match" : "No Match";       
        auto t = TOC_MS(runVerif); 
        writeProcTime(RUNTIMEFILE, t);
        cout << resultStr << endl; 
        cout << "Verification Time : " << t << "ms" << endl; 


    }





    return 0;
}