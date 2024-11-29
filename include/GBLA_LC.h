#ifndef GBLA_LC_H_INCLUDED
#define GBLA_LC_H_INCLUDED

#include "Code.h"
#include "Multinomial.h"

#include <fstream>

using namespace GaloisField;

namespace GBLA_LC
{
    void Codewords(Code &c, wrdout_format wf, std::ostream &out, bool std_out);
    Matrix ReadGFMatrixFromStr(const std::string &str, int rows, int cols);

    //Read a code with parameters stored in a file
    //input: cfile - file that stores codes parameters n,k,H
    //       cnum - number of the code to read (0 -> first code)
    //output: c - code
    bool ReadCodeFromFile(std::ifstream &cfile, std::string &cname, Code &c);


    //Print the Grobner representation of a Code i.e.
    //canonical forms and matphi structure
    //Use the function grobner_representation implemented in the Code class.
    //input: c - Code
    //       wr - word representation (ex. standard, general, etc. See Word.h file)
    //       wm - word metric defined (ex. hamming, rank, etc. See Word.h file)
    //       mord - monomial order (ex. lex, deglex, degrevlex, etc. See Word.h file)
    //       wf - word printed format (ex. as vector, monomial, etc. See Word.h file)
    //
    //Ex. GrobnerRepresentation(c, standard_rep, hamming_m, degrevlex, wf_monomial);
    void GrobnerRepresentation(Code &c, word_representation wr, word_metric wm,
                monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out);


    //Print all coset leaders of a Code
    //Use the function coset_leaders implemented in the Code class.
    //input: assume a standard word representation
    //       for parameters description see the GrobnerRepresentation function
    //
    //Ex. CosetLeaders(cod, hamming_m, degrevlex, wf_monomial_sv);
    void CosetLeaders(Code &c, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out);


    //Print the set of all leader codewords of a Code
    //Use the function leader_codewords implemented in the Code class.
    //input: see parameters description of CosetLeaders function
    //
    //Ex. LeaderCodewords(cod, hamming_m, degrevlex, wf_vector);
    void LeaderCodewords(Code &c, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out);


    bool FindPermutation(const Code &c1, const Code &c2, std::vector<int> &perm,
                            int level=std::numeric_limits<int>::max());

    //**** Experimental functions *********


    //Test the function coset_leaders implemented in the Code class. For each coset leader set
    //obtained with this function the complete coset is generated.
    //input: see parameters description of CosetLeaders function
    //output: true - if all leaders with coset_leaders function match with the leaders in the whole coset generated
    //        false - if a mismatch is found. (the mismatch cosets leaders are printed)
    //
    //Ex. TestCosetLeaders(cod, rank_m, degrevlex, wf_monomial_sv);
    bool TestCosetLeaders(Code &c, word_metric wm,
                    monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out);


    //Compare the coset leaders set returned by the function coset_leaders implemented in the Code class with
    //the GAP coset leaders output using the function Adv_GBLA_CLref_crits3_wt01 implemented in file
    //18022015-Nuevas-IE.txt.
    //input: c - Code
    //       gap_fname - file name with the GAP output.
    //       assume a standard representation.
    //       for the rest parameters description see the GrobnerRepresentation function
    //
    //Ex. CompareCosetLeadersGAP(cod, "C2-10-4cls.txt", hamming_m, degrevlex, wf_monomial_sv);
    void CompareCosetLeadersGAP(Code &c, const std::string &gap_fname, word_metric wm, monomial_order mord, wrdout_format wf);

    void AutomorphismGroup(const Code &cod, const Code &pcod, int level, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out);
    void PermutationAnalysis(const Code &cod, const Code &pcod, int level, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out, bool detailed);

    void LeaderCodewordsByLevels(const Code &c, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out);

    void PermuteCode(Code &c, std::vector<int> &perm,
                     std::ostream &out, bool std_out);


    //Output utility functions
    void ClearScreen();
    void SystemPause();
}


#endif // GBLA_LC_H_INCLUDED
