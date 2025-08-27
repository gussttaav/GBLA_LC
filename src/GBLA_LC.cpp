#include "GBLA_LC.h"
#include <time.h>
#include <sstream>
#include <stdlib.h>

#ifdef _WIN32
    // Windows-specific code
    #define CLEAR_CMD "cls"
    #define PAUSE_CMD "PAUSE"
#else
    // Unix/Linux/macOS-specific code
    #define CLEAR_CMD "clear"
    #define PAUSE_CMD "echo 'Press Enter to continue...'; read line"
#endif

using namespace GaloisField;

namespace GBLA_LC
{
    //**** utils

    void ClearScreen(){
        if(system(CLEAR_CMD) != 0)
            std::cout << "\n\n";
    }

    void SystemPause(){
        if(system(PAUSE_CMD) != 0)
            std::cout << "\n\n";
    }

    void ShowCodeParameters(const Code &c, std::ostream &out)
    {
        out << "[" << c.get_length() << ", " << c.get_dimension()
                  << "] Linear code over GF(" << GaloisField::q << ")\n"
                  << "Generator Matrix: \n" << c.get_gen_mat()
                  << "\nCheck Matrix: \n" << c.get_check_mat() << "\n\n";
    }


    void ShowGBLAComputingParameters(word_metric wm,
                monomial_order mo, std::ostream &out)
    {
        if(wm == hamming_m)
            out << "hamming";
        else out << "rank";
        out << " word metric and a ";
        switch(mo){
            case lex: out << "lexicographical"; break;
            case deglex: out << "degree lexicographical"; break;
            case degrevlex: out << "degree reverse lexicographical"; break;
        }
        out << " order.";
    }


    template<typename T>
    void PrintVector(const std::vector<T> &vec, std::ostream &out)
    {
        if(vec.empty())
        {
            out << "[]";
            return;
        }

        typename std::vector<T>::const_iterator i = vec.begin();
        out << "[" << *i;
        for(i++; i != vec.end(); i++)
            out << ", " << *i;
        out << "]";
    }

    template<typename T>
    void PrintVector(const std::set<T> &vec, std::ostream &out)
    {
        if(vec.empty())
        {
            out << "[]";
            return;
        }

        typename std::set<T>::const_iterator i = vec.begin();
        out << "[" << *i;
        for(i++; i != vec.end(); i++)
            out << ", " << *i;
        out << "]";
    }

    template<typename T>
    void PrintVector(const std::list<T> &vec, std::ostream &out)
    {
        if(vec.empty())
        {
            out << "[]";
            return;
        }

        typename std::list<T>::const_iterator i = vec.begin();
        out << "[" << *i;
        for(i++; i != vec.end(); i++)
            out << ", " << *i;
        out << "]";
    }

    template<typename T>
    void PrintMatrix(const std::list<std::set<T> > &mat, std::ostream &out)
    {
        typename std::list<std::set<T> >::const_iterator i = mat.begin();
        out << "[ ";
        PrintVector(*i, out);
        for(i++; i != mat.end(); i++)
        {
            out << ", ";
            PrintVector(*i, out);
        }
        out << " ]";
    }

    template<typename T>
    void PrintMatrix(const std::list<std::list<T> > &mat, std::ostream &out)
    {
        typename std::list<std::list<T> >::const_iterator i = mat.begin();
        out << "[ ";
        PrintVector(*i, out);
        for(i++; i != mat.end(); i++)
        {
            out << ", ";
            PrintVector(*i, out);
        }
        out << " ]";
    }

    template<typename T>
    void PrintMatrix(const std::list<std::vector<T> > &mat, std::ostream &out)
    {
        typename std::list<std::vector<T> >::const_iterator i = mat.begin();
        out << "[ ";
        PrintVector(*i, out);
        for(i++; i != mat.end(); i++)
        {
            out << ", ";
            PrintVector(*i, out);
        }
        out << " ]";
    }

    Matrix ReadGFMatrixFromStr(const std::string &str, int rows, int cols)
    {
        int sti = 1, ste = 0;
        std::string st_num;
        Matrix h(rows, cols);

        for(int i = 0; i < rows; i++)
        {
            while(str[sti] != '[') sti++;
            ste = sti;
            for(int j = 0; j < cols; j++)
            {
                sti = ste+1;
                ste = sti;
                while(str[ste] != ',' && str[ste] != ']') ste++;
                st_num = str.substr(sti, ste-sti);
                h(i,j) = GFElem(st_num);
            }
            sti = ste+1;
        }

        return h;
    }


    void Codewords(Code &c, wrdout_format wf, std::ostream &out, bool std_out)
    {
        std::set<Word> cwords;
        wrdout_format ftmp = Word::get_wrdout_format();
        Word::set_wrdout_format(wf);

        ClearScreen();
        ShowCodeParameters(c, std::cout);
        std::cout << "\nGenerating codewords. Wait...";
        clock_t exe_time = clock();

        c.gen_codewords(cwords);
        exe_time = clock() - exe_time;
        ClearScreen();
        ShowCodeParameters(c, out);
        out << "\n*** CODEWORDS ***\n\n"
                  << "Execution time: " << std::setprecision(7) << exe_time << " ("
                  << ((double)exe_time/CLOCKS_PER_SEC) << " seconds)\n\n";
        PrintVector(cwords, out);
        out << "\n\n";

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();

        Word::set_wrdout_format(ftmp); //restore word output format*/
    }

    //Print the Grobner representation of a Code i.e.
    //canonical forms and matphi structure
    void GrobnerRepresentation(Code &c, word_representation wr, word_metric wm,
                monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out)
    {
        std::vector<Word> X, N;
        std::set<matphi_elem> matphi;
        wrdout_format ftmp = Word::get_wrdout_format();
        Word::set_wrdout_format(wf);

        ClearScreen();
        ShowCodeParameters(c, std::cout);
        std::cout << "\nComputing the Grobner Representation using the ";
        ShowGBLAComputingParameters(wm, mord, std::cout);
        std::cout << " Wait...";
        clock_t exe_time = clock();

        c.grobner_representation(N, matphi, wr, wm, mord);
        exe_time = clock() - exe_time;
        ClearScreen();
        ShowCodeParameters(c, out);
        out << "\n*** GROBNER REPRESENTATION ***\nComputed with ";
        ShowGBLAComputingParameters(wm, mord, out);
        out << "\n\n" << "Execution time: " << std::setprecision(7) << exe_time
                  << " ("  << ((double)exe_time/CLOCKS_PER_SEC) << " seconds)\n\n";

        out << "*** GBLA variables\n\n";
        c.generate_variables(X, wr);
        out << X.size() << " elements\n";
        PrintVector(X, out);

        out << "\n\n*** CANONICAL FORMS\n\n"
                  << N.size() << " elements\n";
        PrintVector(N, out);

        out << "\n\n*** MATPHI STRUCTURE\n\n"
                  << matphi.size() << " elements\n";
        std::set<matphi_elem>::iterator i = matphi.begin();
        out << "[ [" << N[(*i).dom_N] << ", "
                  << X[(*i).dom_X] << ", " << N[(*i).img_N] << "]";
        for(i++; i != matphi.end(); i++)
            out << ", [" << N[(*i).dom_N] << ", "
                  << X[(*i).dom_X] << ", " << N[(*i).img_N] << "]";
        out << " ]\n\n";

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();

        Word::set_wrdout_format(ftmp); //restore word output format*/
    }



    //Print all coset leaders of a Code
    //Use the function coset_leaders implemented in the Code class.
    void CosetLeaders(Code &c, word_metric wm, monomial_order mord,
                wrdout_format wf, std::ostream &out, bool std_out)
    {
        std::vector<Word> X;
        std::list<std::set<Word> > CL;
        wrdout_format ftmp = Word::get_wrdout_format();
        Word::set_wrdout_format(wf);

        ClearScreen();
        ShowCodeParameters(c, std::cout);
        std::cout << "\nComputing the cosets leaders using the ";
        ShowGBLAComputingParameters(wm, mord, std::cout);
        std::cout << " Wait...";
        clock_t exe_time = clock();

        c.coset_leaders(CL, wm, mord);
        exe_time = clock() - exe_time;
        ClearScreen();
        ShowCodeParameters(c, out);
        out << "\n*** COSETS LEADERS USING GBLA ALGORITHM ***\nComputed with ";
        ShowGBLAComputingParameters(wm, mord, out);
        out << "\n\n" << "Execution time: " << std::setprecision(7) << exe_time
                  << " ("  << ((double)exe_time/CLOCKS_PER_SEC) << " seconds)\n\n";

        out << "*** GBLA variables\n\n";
        c.generate_variables(X, standard_rep);
        out << X.size() << " elements\n";
        PrintVector(X, out);

        out << "\n\n*** COSETS LEADERS\n"
                  << CL.size() << " cosets\n\n";
        PrintMatrix(CL, out);
        out << "\n\n";

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();

        Word::set_wrdout_format(ftmp); //restore word output format*/
    }



    //Print the set of all leader codewords of a Code
    //Use the function leader_codewords implemented in the Code class.
    void LeaderCodewords(Code &c, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out)
    {
        std::vector<Word> X;
        std::set<Word> lcwords;
        wrdout_format ftmp = Word::get_wrdout_format();
        Word::set_wrdout_format(wf);

        ClearScreen();
        ShowCodeParameters(c, std::cout);
        std::cout << "\nComputing the leader codewords according "
                  << "to ZNTS_AAECC Def. 5 p. 6 using the ";
        ShowGBLAComputingParameters(wm, mord, std::cout);
        std::cout << " Wait...";
        clock_t exe_time = clock();

        c.leader_codewords(lcwords, wm, mord);
        exe_time = clock() - exe_time;
        ClearScreen();
        ShowCodeParameters(c, out);
        out << "\n*** LEADER CODEWORDS USING GBLA ALGORITHM "
                  << "(ZNTS_AAECC Def. 5 p. 6) ***\nComputed with ";
        ShowGBLAComputingParameters(wm, mord, out);
        out << "\n\nExecution time: " << std::setprecision(7) << exe_time << " ("
                  << ((double)exe_time/CLOCKS_PER_SEC) << " seconds)\n\n";

        out << "*** GBLA variables\n\n";
        c.generate_variables(X, standard_rep);
        out << X.size() << " elements\n";
        PrintVector(X, out);

        out << "\n\n*** LEADER CODEWORDS\n"
                  << lcwords.size() << " elements \n\n";
        PrintVector(lcwords, out);
        out << "\n\n";

        /*double gap_exe_time = 138891;
        std::cout << "Execution time: " << gap_exe_time
                  << " ("  << ((double)gap_exe_time/CLOCKS_PER_SEC) << " seconds)\n\n";*/
        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();

        Word::set_wrdout_format(ftmp); //restore word output format*/
    }



    //Error correct capability using the coset leaders set
    //Find the first coset with more than one leader.
    //error-correcting-capability = coset-weight - 1.
    int ErrorCorrectCapability(std::list<std::set<Word> > CL)
    {
        std::list<std::set<Word> >::iterator itl;

        for(itl = CL.begin(); itl != CL.end(); itl++)
        {
            if((*itl).size() != 1)
                return (*((*itl).begin())).get_weight()-1;
        }

        itl = CL.begin();
        itl++;
        return (*((*itl).begin())).get_weight()-1;
    }


    //Test the function coset_leaders implemented in the Code class. For each coset leader set
    //obtained with this function the complete coset is generated.
    bool TestCosetLeaders(Code &c, word_metric wm,
                monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out)
    {
        int cnum = 1;
        bool match = true;
        std::set<Word> coset;
        std::list<std::set<Word> > CL;
        std::set<Word>::iterator its, itc;
        std::list<std::set<Word> >::iterator itl;

        ClearScreen();
        if(!std_out)
        {
            ShowCodeParameters(c, std::cout);
            std::cout << "\n*** TESTING GBLA ALGORITHM FOR COSETS LEADERS***\nComputed with ";
            ShowGBLAComputingParameters(wm, mord, std::cout);
            std::cout << "\n(All cosets are generated exhaustively. In each coset words are "
                  << "ordered by weight and the leaders must be found in the coset leaders set "
                  << "computed using GBLA.)\n\n";
        }

        ShowCodeParameters(c, out);
        out << "\n*** TESTING GBLA ALGORITHM FOR COSETS LEADERS***\nComputed with ";
        ShowGBLAComputingParameters(wm, mord, out);
        out << "\n(All cosets are generated exhaustively. In each coset words are "
                  << "ordered by weight and the leaders must be found in the coset leaders set "
                  << "computed using GBLA.)\n\n";

        c.coset_leaders(CL, wm, mord);

        wrdout_format ftmp = Word::get_wrdout_format();
        Word::set_wrdout_format(wf);

        out << "Error-correcting capability = "
                  << ErrorCorrectCapability(CL) << "\n\n";

        for(itl = CL.begin(); itl != CL.end(); itl++, cnum++)
        {
          its = (*itl).begin();
          c.gen_coset(*its, coset);
          itc = coset.begin();
          std::advance(itc, (*itl).size());

          if((*itc).get_weight() == (*its).get_weight())
          {  match = false;  break;  }
          else{
            out << "COSET " << cnum << "  (weight="
                      << (*coset.begin()).get_weight() << ")\n";
            PrintVector(coset, out);
            out << "\n\n";
          }
        }

        if(!match)
        {
            out << "*** A COSET LEADERS MISMATCH WERE FOUND\n\n"
                      << "COSET LEADER\n";
            PrintVector(*itl, out);
            out << "\n\nCOMPLETE COSET\n";
            PrintVector(coset, out);
            out << "\n\n";
        }
        else {
            out << "*** COSETS LEADERS BY GBLA\n\n";
            PrintMatrix(CL, out);
            out << "\n\n*** ALL COSET LEADERS WERE FOUND ***\n\n";
        }

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();

        Word::set_wrdout_format(ftmp); //restore word output format*/
        return match;
    }



    bool ReadCodeFromFile(std::ifstream &cfile, std::string &cname, Code &c)
    {
        if(!cfile.is_open()) return false;

        char nbuffer[15];
        std::size_t found;
        bool code_found = false;
        int n, nend, forder = 2, cn = 1, ck = 1;
        std::string line, nline, st_num, str_alpha;

        GFElem alpha_elem;

        while(std::getline(cfile, line))
        {
            n = line.length();
            found = line.find("-->");

            //--> GF[order], [alpha]
            if(found != std::string::npos)
            {
                //Read Galois Field order
                for(nend = 6; isdigit(line[nend]); nend++);
                st_num = line.substr(6, nend-6);
                strcpy(nbuffer, st_num.c_str());
                forder = atoi(nbuffer);

                //Read primitive element
                if(nend < n && line[nend] == ',')
                {
                    str_alpha = line.substr(nend+1);
                    alpha_elem = GFElem(str_alpha);
                }
                else str_alpha = "1";

                continue;
            }

            if(line.find(cname)!= std::string::npos)
                code_found = true;

            found = line.find("n=");
            if(found != std::string::npos)
            {
                //read code length
                for(nend = found+2; isdigit(line[nend]); nend++);
                st_num = line.substr(found+2, nend-(found+2));
                strcpy(nbuffer, st_num.c_str());
                cn = atoi(nbuffer);
            }

            found = line.find("k=");
            if(found != std::string::npos)
            {
                //read code dimension
                for(nend = found+2; isdigit(line[nend]); nend++);
                st_num = line.substr(found+2, nend-(found+2));
                strcpy(nbuffer, st_num.c_str());
                ck = atoi(nbuffer);
                continue;
            }

            found = line.find("[");
            if(found != std::string::npos && code_found)
            {
                GF::set_order(forder, str_alpha);
                Matrix h = ReadGFMatrixFromStr(line, cn-ck, cn);

                std::getline(cfile, nline);
                if(nline.length() > 1)
                {
                    Matrix g = ReadGFMatrixFromStr(nline, ck, cn);
                    c = Code(g, h);
                }
                else c = Code(h);

                return true;
            }
        }

        return false;
    }


    //begin **** Functions to read GAP output
    #define LEFT_BRK    1
    #define RIGHT_BRK   2
    #define COMMA_SEP   3
    #define VAR_NAME    4
    #define VAR_NUM     5
    #define VAR_PROD    6
    #define ERROR       7
    #define END_STR     8
    #define VAR_EXP     9

    int next_token(const std::string &str, std::string::size_type &pos, int &num)
    {
        char nbuffer[256];
        std::string st_num;
        std::string::size_type ipos;

        while(str[pos] == ' ') pos++;
        if(pos >= str.size()) return END_STR;

        if(isdigit(str[pos]))
        {
            ipos = pos;
            while(isdigit(str[pos])) pos++;
            st_num = str.substr(ipos, pos-ipos);
            strcpy(nbuffer, st_num.c_str());
            num = atoi(nbuffer);
            return VAR_NUM;
        }

        int token;
        switch(str[pos])
        {
            case '[': token = LEFT_BRK;
                break;
            case ']': token = RIGHT_BRK;
                break;
            case ',': token = COMMA_SEP;
                break;
            case '*': token = VAR_PROD;
                break;
            case 'x': if(pos+1 < str.size() && str[pos+1] == '_')
                      {
                          token = VAR_NAME;
                          pos++;
                      }
                      else token = ERROR;
                break;
            case '^': 
                    for(ipos = pos+1; ipos < str.size() && isdigit(str[ipos]); ipos++);
                    if(pos+1 != ipos)
                    {
                        token = VAR_EXP;
                        st_num = str.substr(pos+1, ipos-pos-1);
                        strcpy(nbuffer, st_num.c_str());
                        num = atoi(nbuffer);
                        pos = ipos-1;
                    }
                    else token = ERROR;
                break;
                
            default:
                    token = ERROR;
                break;
        }
        pos++;
        return token;
    }

    bool ReadCosetLeadersFromFileGAP(std::ifstream &clsfile, std::list<std::set<Word> > &leaders, double &gap_exe_time)
    {
        Word nw;
        int token, num, vnum, vexp;
        leaders.clear();
        char nbuffer[256];
        std::set<Word> nset;
        std::string cls, line;
        std::string::size_type pos = 0;

        //Read execution time
        std::getline(clsfile, line);
        for(std::string::size_type i = 0; i < line.size(); i++)
            if(!isdigit(line[i]))
               return false;
        strcpy(nbuffer, line.c_str());
        gap_exe_time = atof(nbuffer);

        //read coset leaders as string type
        while(!clsfile.eof())
        {
            std::getline(clsfile, line);
            cls += line;
        }

        //read empty word at first
        //seq [1],[
        for(int i = 0; i < 5; i++)
           token = next_token(cls, pos, num);
        nset.insert(Word());
        leaders.push_back(nset);

        while(token != END_STR && token != ERROR)
        {
            //read new Coset Leader
            if(token == LEFT_BRK)
            {
                nset.clear();
                token = next_token(cls, pos, num);

                while(token != END_STR && token != COMMA_SEP)
                {
                    nw = Word();
                    while(token != COMMA_SEP && token != RIGHT_BRK)
                    {
                        vnum = -1;
                        vexp = -1;

                        if(token == VAR_NUM)
                            vnum = num;

                        if(token == END_STR || token == ERROR)
                           return false;

                        token = next_token(cls, pos, num);

                        if(token == VAR_EXP)
                            vexp = num;

                        if(vnum != -1)
                        {
                            if(vexp == -1)
                               nw.set_monomial_var_exponent(vnum-1, 1);
                            else nw.set_monomial_var_exponent(vnum-1, vexp);
                        }
                    }

                    nset.insert(nw);
                    token = next_token(cls, pos, num);
                }

                leaders.push_back(nset);
            }

            token = next_token(cls, pos, num);
        }

        if(token == ERROR)
            return false;
        return true;
    }

    //end **** Functions to read GAP output


    //Compare the coset leaders set returned by the function coset_leaders implemented in the Code class with
    //the GAP coset leaders output using the function Adv_GBLA_CLref_crits3_wt01 implemented in file
    //18022015-Nuevas-IE.txt.
    void CompareCosetLeadersGAP(Code &c, const std::string &gap_fname, word_metric wm,
                    monomial_order mord, wrdout_format wf = wf_monomial)
    {
        double gap_exe_time;
        bool found_coset, found_word;
        std::set<Word>::iterator gits, cits;
        std::list<std::set<Word> >::iterator gitl, citl;
        std::list<std::set<Word> > gap_leaders, cleaders;
        std::ifstream clsfile(gap_fname.c_str());
        wrdout_format ftmp = Word::get_wrdout_format();
        Word::set_wrdout_format(wf);

        ShowCodeParameters(c, std::cout);
        std::cout << "Building Coset Leaders using GBLA ...\n";
        clock_t exe_time = clock();
        c.coset_leaders(cleaders, wm, mord);
        exe_time = clock() - exe_time;
        std::cout << "Coset Leaders completed. Execution time: " << std::setprecision(7)
                  << exe_time << " ("  << ((double)exe_time/CLOCKS_PER_SEC) << " seconds)\n\n"
                  << "Reading GAP output for Coset Leaders using GBLA ...\n";
        ReadCosetLeadersFromFileGAP(clsfile, gap_leaders, gap_exe_time);
        std::cout << "GAP Coset Leaders completed. Execution time: " << gap_exe_time
                  << " (" << gap_exe_time/CLOCKS_PER_SEC << " seconds)\n\n";

        for(gitl = gap_leaders.begin(); gitl != gap_leaders.end(); gitl++)
        {
            found_coset = false;
            found_word = true;
            gits = (*gitl).begin();

            //find GAP coset leader in cleader
            for(citl = cleaders.begin(); citl != cleaders.end(); citl++ )
            {
                cits = (*citl).find(*gits);
                if(cits != (*citl).end())
                {
                    found_coset = true;
                    break;
                }
            }

            if(found_coset)
            {
                //Found coset leader
                //Now compare words of coset leaders
                for(gits = (*gitl).begin(); gits != (*gitl).end(); gits++)
                {
                    if((*citl).find(*gits) == (*citl).end())
                    {
                        found_word = false;
                        break;
                    }
                }

                if(!found_word)
                    break;
            }
            else break;
        }

        if(!found_coset || !found_word)
        {
            if(!found_coset)
                std::cout << "ERROR !!! GAP COSET LEADER NOT FOUND\n\n";
            else {
                std::cout << "ERROR !!! COSET LEADERS MISMATCH\n\n"
                          << "COSET LEADER USING C++ IMPLEMENTATION\n";
                PrintVector(*citl, std::cout);
                std::cout << "\n\nCOSET LEADER USING GAP IMPLEMENTATION\n[";
            }

            PrintVector(*gitl, std::cout);
            std::cout << "\n\n";
        }
        else std::cout << "*** ALL COSETS MATCH PERFECTLY *** \n\n";

        std::cout << "COSET LEADERS USING C++ IMPLEMENTATION\n"
                  << cleaders.size() << " cosets \n";
        PrintMatrix(cleaders, std::cout);
        std::cout << "\n\nCOSET LEADERS USING GAP IMPLEMENTATION\n"
                  << gap_leaders.size() << " cosets \n";
        PrintMatrix(gap_leaders, std::cout);

        Word::set_wrdout_format(ftmp); //restore word output format
    }




    //**** Code equivalence functions ***
    //

    typedef std::vector<int> Interaction;
    typedef std::vector<std::pair<int, Interaction> > TagItem;
    typedef std::vector<TagItem> Tag;
    typedef std::vector<Tag> Label;

    typedef std::pair<std::pair<Tag, int>, std::list<int> > SigItem;
    typedef std::list<SigItem> Signature;

    void PartitionWeight(const std::set<Word> &LC, std::list<std::set<Word> > &part, int level)
    {
        part.clear();
        if(LC.empty()) return;

        std::set<Word> part_elem;
        std::set<Word>::const_iterator it = LC.begin();
        int part_weight, lc = 0;
        part_weight = it->get_weight() == 0 ? (++it)->get_weight() : it->get_weight();

        while(it != LC.end() && lc < level)
        {
            if(it->get_weight() == part_weight)
                part_elem.insert(*it);
            else{
                part.push_back(part_elem);
                part_elem.clear();
                part_elem.insert(*it);
                part_weight = it->get_weight();
                lc++;
            }
            it++;
        }

        if(it == LC.end() && lc < level && !part_elem.empty())
            part.push_back(part_elem);
    }

    void print_tag_item(const TagItem &t, int level)
    {
        if(t.empty())
        {
            std::cout << "0";
            return;
        }

        TagItem::size_type i = 0;
        Interaction::size_type j = 0;
        if(t[0].first != 1)
            std::cout << t[0].first << "*";
        while(t[0].second[j] == 0) j++;
        if(j == t[0].second.size()-1)
        {
            std::cout << "z" << Word::get_length()+1;
            if(t[0].second[j] > 1) std::cout << "^" << t[0].second[j];
        }
        else{
            std::cout << "z" << j+1;
            if(t[0].second[j] > 1) std::cout << "^" << t[0].second[j];
            for(j++; j < t[0].second.size()-1; j++)
                if(t[0].second[j] != 0)
                {
                    std::cout << "*z" << j+1;
                    if(t[0].second[j] > 1) std::cout << "^" << t[0].second[j];
                }

            std::cout << "*z" << Word::get_length()+1;
            if(t[0].second[j] > 1) std::cout << "^" << t[0].second[j];
        }

        for(i++; i < t.size(); i++)
        {
            j = 0;
            std::cout << " + ";
            if(t[i].first != 1)
                std::cout << t[i].first << "*";
            while(t[i].second[j] == 0) j++;
            if(j == t[i].second.size()-1)
            {
                std::cout << "z" << Word::get_length()+1;
                if(t[i].second[j] > 1) std::cout << "^" << t[i].second[j];
            }
            else{
                std::cout << "z" << j+1;
                if(t[i].second[j] > 1) std::cout << "^" << t[i].second[j];
                for(j++; j < t[i].second.size()-1; j++)
                    if(t[i].second[j] != 0)
                    {
                        std::cout << "*z" << j+1;
                        if(t[i].second[j] > 1) std::cout << "^" << t[i].second[j];
                    }
                std::cout << "*z" << Word::get_length()+1;
                if(t[i].second[j] > 1) std::cout << "^" << t[i].second[j];
            }
        }
    }

    void print_label(const Label &lb, int level)
    {
        Label::const_iterator i = lb.begin();
        Tag::const_iterator j = i->begin();
        std::cout << "[ [";
        print_tag_item(*j, level);
        for(j++; j != i->end(); j++)
        {
            std::cout << ", ";
            print_tag_item(*j, level);
        }
        std::cout << " ]";

        for(i++; i != lb.end(); i++)
        {
            Tag::const_iterator j = i->begin();
            std::cout << ", [ ";
            print_tag_item(*j, level);
            for(j++; j != i->end(); j++)
            {
                std::cout << ", ";
                print_tag_item(*j, level);
            }
            std::cout << " ]";
        }
        std::cout << " ]";
    }


    void LabelMonomialsByLevel(const std::set<Word> &l, const std::list<int> &S, Label &lb, int level)
    {
        int z;
        std::list<int>::iterator i;
        std::list<int>::const_iterator j;
        std::set<Word>::const_iterator li;
        std::list<int> sp, fsp, isp, u = S;
        u.sort();
        Interaction sg_inter(S.size()+1, 0), inter;

        for(li = l.begin(); li != l.end(); li++)
        {
            fsp.clear();
            isp.clear();
            inter = sg_inter;

            sp = li->support();
            std::set_difference(sp.begin(), sp.end(),
                    u.begin(), u.end(), std::inserter(fsp, fsp.begin()));

            //callTerm function
            std::set_intersection(sp.begin(), sp.end(),
                    u.begin(), u.end(), std::inserter(isp, isp.begin()));
            for(i = isp.begin(); i != isp.end(); i++)
                for(z = S.size()-1, j = S.begin(); j != S.end(); j++, z--)
                    if(*i == *j){ inter[z] = (int)(*li)[*i]; break; }

            //add interaction
            bool found_interaction;
            for(i = fsp.begin(); i != fsp.end(); i++)
            {
                found_interaction = false;
                inter[S.size()] = (int)(*li)[*i];
                for(TagItem::size_type ti = 0; ti < lb[*i][level].size() && !found_interaction; ti++)
                    if(lb[*i][level][ti].second == inter)
                    {
                        lb[*i][level][ti].first++;
                        found_interaction = true;
                    }
                if(!found_interaction)
                    lb[*i][level].push_back(std::make_pair(1, inter));
            }
        }
    }


    void LabelMonomials(const std::list<std::set<Word> > &P, const std::list<int> &S, Label &lb)
    {
        int level_i;
        std::list<std::set<Word> >::const_iterator Pi;
        lb = Label(Word::get_length(), Tag(P.size()));

        for(level_i = 0, Pi = P.begin(); Pi != P.end(); Pi++, level_i++)
            LabelMonomialsByLevel(*Pi, S, lb, level_i);
        //print_label(lb, P.size());
        //std::cout << "\n\n";
    }

    bool EqualTags(const Tag &t1, const Tag &t2)
    {
        bool item_found;
        for(Tag::size_type i = 0; i < t1.size(); i++)
        {
            if(t1[i].size() != t2[i].size())
                return false;

            for(TagItem::size_type j = 0; j < t1[i].size(); j++)
            {
                item_found = false;
                for(TagItem::size_type k = 0; k < t1[i].size(); k++)
                    if(t1[i][j] == t2[i][k])
                    {   item_found = true;  break;  }

                if(!item_found)
                    return false;
            }
        }

        return true;
    }

    bool EmptyTag(const Tag &t)
    {
        for(Tag::const_iterator i = t.begin(); i != t.end(); i++)
            if(!i->empty()) return false;
        return true;
    }

    bool CmpSignature(const SigItem &s1, const SigItem &s2)
    {
        if(s1.first.second != s2.first.second)
            return s1.first.second < s2.first.second;

        int wl, cnum, cmp = -1;
        wl = Word::get_length();
        Word::set_length(wl+1);

        Term et, t1, t2;
        Multinomial ep, p1, p2;
        TagItem::const_iterator j1, j2;
        Interaction::const_iterator k1, k2;
        Tag::const_iterator i1 = s1.first.first.begin(), i2 = s2.first.first.begin();

        for(; i1 != s1.first.first.end() && cmp == -1; i1++, i2++)
        {
            if(i1->empty() || i2->empty())
            {
                if(i1->empty() && !i2->empty())
                    cmp = 1; //true
                else if(!i1->empty() && i2->empty())
                    cmp = 0; //false
            }
            else{
                //Construct multinomial associated
                p1 = ep;  p2 = ep;
                for(j1 = i1->begin(); j1 != i1->end(); j1++)
                {
                    t1 = et;
                    k1 = j1->second.begin();
                    for(int i = 0; k1 != j1->second.end(); i++, k1++)
                        t1[i] = *k1;
                    t1.set_coefficient(j1->first);
                    p1 += t1;
                }
                for(j2 = i2->begin(); j2 != i2->end(); j2++)
                {
                    t2 = et;
                    k2 = j2->second.begin();
                    for(int i = 0; k2 != j2->second.end(); i++, k2++)
                        t2[i] = *k2;
                    t2.set_coefficient(j2->first);
                    p2 += t2;
                }

                //Compare multinomial
                if(p1.similar(p2))
                {
                    //compare by coefficients
                    cnum = p1.terms_num();
                    for(int i = 0; i < cnum && cmp == -1; i++)
                        if(p1[i].get_coefficient() != p2[i].get_coefficient())
                            cmp = p1[i].get_coefficient() < p2[i].get_coefficient() ? 1 : 0;
                }
                else cmp = p1 < p2 ? 1 : 0;
            }
        }

        Word::set_length(wl);
        return cmp == 1 ? true : false;
    }


    bool EqualSignature(const Signature &s1, const Signature &s2)
    {
        if(s1.size() != s2.size())
            return false;

        Signature::const_iterator i, j;
        for(i = s1.begin(), j = s2.begin(); i != s1.end(); i++, j++)
            if(i->second.size() != j->second.size())
                return false;

        return true;
    }

    void ListSignature(const Label &lb, const std::list<int> &S, Signature &sg)
    {
        sg.clear();
        int pos = 0;
        bool tag_found;
        Signature::iterator si;

        for(Label::const_iterator i = lb.begin(); i != lb.end(); i++, pos++)
        {
            if(std::find(S.begin(), S.end(), pos) == S.end())
            {
                //Find Tag in the Label list
                tag_found = false;
                for(si = sg.begin(); si != sg.end(); si++)
                {
                    if(EqualTags(*i, si->first.first))
                    {
                        si->first.second++;
                        tag_found = true;
                        si->second.push_back(pos);
                        break;
                    }
                }
                if(!tag_found/* && !EmptyTag(*i)*/)
                   sg.push_back(std::make_pair(std::make_pair(*i, 1), std::list<int>(1, pos)));
            }
        }

        sg.sort(CmpSignature);

        //for debugging purpose only
        /*std::cout << "[";
        for(si = sg.begin(); si != sg.end(); si++)
        {
            Tag::iterator ti = si->first.first.begin();
            std::cout << "[[";
            print_tag_item(*ti, 2);
            for(ti++; ti != si->first.first.end(); ti++)
            {
                std::cout << ", ";
                print_tag_item(*ti, 2);
            }
            std::cout << " ], " << si->first.second << "], [";

            std::list<int>::iterator li = si->second.begin();
            std::cout << *li;
            for(li++; li != si->second.end(); li++)
                std::cout << ", " << *li;
            std::cout << "]], ";
        }

        std::cout << "\n\n";*/
    }


    bool IsValidPermutation(const std::list<std::set<Word> > &ll1,
                    const std::list<std::set<Word> > &ll2, const std::vector<int> &perm)
    {
        Word w;
        std::set<Word>::const_iterator j;
        std::list<std::set<Word> >::const_iterator i1, i2;

        for(i1 = ll1.begin(), i2 = ll2.begin(); i1 != ll1.end(); i1++, i2++)
        {
            for(j = i1->begin(); j != i1->end(); j++)
            {
                w = *j;
                w.permute(perm);
                if(i2->find(w) == i2->end())
                    return false;
            }
        }

        return true;
    }

    void print_signature(const Signature &sg, std::ostream &out)
    {
        Signature::const_iterator i = sg.begin();
        out << "[";
        PrintVector(i->second, out);
        for(i++; i != sg.end(); i++)
        {
            out << ", ";
            PrintVector(i->second, out);
        }
        out << "]";
    }

    void PrintPermutationAsCycles(const std::vector<int> &perm, std::ostream &out)
    {
        int len = perm.size(), j;
        std::vector<bool> visited(len, false);
        std::string pcycles;
        std::stringstream cycle;

        for(int i = 0; i < len; i++)
        {
            if(!visited[i] && perm[i] != i)
            {
                cycle.str("");
                cycle << "(" << i;
                j = perm[i];

                while(j != i)
                {
                    cycle << "," << j;
                    visited[j] = true;
                    j = perm[j];
                }

                cycle << ")";
                pcycles += cycle.str();
            }
            visited[i] = true;
        }

        if(!pcycles.empty())
            out << pcycles;
        else out << "()";
    }

    bool is_identity_permutation(const std::vector<int> &perm)
    {
        int len = perm.size();
        for(int i = 0; i < len; i++)
        {
            if(perm[i] != i)
                return false;
        }
        return true;
    }

    bool FindPermutation(const std::list<std::set<Word> > &ll1, const std::list<std::set<Word> > &ll2,
                std::list<int> &S, std::list<int> &P, std::list<std::vector<int> > &pList, int tl,
                clock_t bg, std::ostream &out, bool print_details)
    {
        if(P.size() == (unsigned)Word::get_length())
        {
            std::list<int>::const_iterator Si, Pi;
            std::vector<int> perm(Word::get_length());
            for(Si = S.begin(), Pi = P.begin(); Si != S.end(); Si++, Pi++)
                perm[*Si] = *Pi;
            pList.push_back(perm);
            if((pList.size() == 1 && !is_identity_permutation(*pList.begin()))
                || pList.size() == 2)
            {
                std::list<std::vector<int> >::iterator pl = pList.begin();
                if(pList.size() == 2) pl++;

                out << "First permutation found: ";
                PrintPermutationAsCycles(*pl, out);
                out << "\nTime: " << std::setprecision(7)
                    << (double)(clock()-bg)/CLOCKS_PER_SEC << " seconds.\n\n";
                return true;
            }

            return false;
        }

        Label f, g;
        Signature lf, lg;

        LabelMonomials(ll1, S, f);
        LabelMonomials(ll2, P, g);
        ListSignature(f, S, lf);
        ListSignature(g, P, lg);

        if(EqualSignature(lf, lg))
        {
            std::list<int>::reverse_iterator j;
            Signature::iterator lfi = lf.begin(), lgi = lg.begin();

            if(lf.back().second.size() == 1)
            {
                int num_ins = 0;
                if(print_details)
                {
                    out << "Tree level: " << tl << "\nlf := ";
                    print_signature(lf, out);
                    out << "\nlg := ";
                    print_signature(lg, out);
                    out << "\n";
                }

                for(; lfi != lf.end(); lfi++, lgi++, num_ins++)
                {
                    S.push_front(lfi->second.front());
                    P.push_front(lgi->second.front());

                    if(print_details)
                       out << S.front() << " -> " << P.front() << "\n";
                }

                if(print_details)
                {
                    out << "S := ";
                    PrintVector(S, out);
                    out << "\nP := ";
                    PrintVector(P, out);
                    out << "\n\n";
                }

                if(!FindPermutation(ll1, ll2, S, P, pList, tl+1, bg, out, print_details))
                {
                    for(int i = 0; i < num_ins; i++)
                    {  S.pop_front();  P.pop_front(); }
                }
                else return true;
            }
            else {

                S.push_front(lfi->second.back());

                for(j = lgi->second.rbegin(); j != lgi->second.rend(); ++j)
                {
                    P.push_front(*j);

                    if(print_details)
                    {
                        out << "Tree level: " << tl << "\nlf := ";
                        print_signature(lf, out);
                        out << "\nlg := ";
                        print_signature(lg, out);
                        out << "\n" << S.front() << " -> " << P.front() << "\nS := ";
                        PrintVector(S, out);
                        out << "\nP := ";
                        PrintVector(P, out);
                        out << "\n\n";
                    }

                    if(!FindPermutation(ll1, ll2, S, P, pList, tl+1, bg, out, print_details))
                        P.pop_front();
                    else return true;
                }

                S.pop_front();
            }
        }

        return false;
    }

    void FindAllPermutations(const std::list<std::set<Word> > &ll1, const std::list<std::set<Word> > &ll2,
                std::list<int> &S, std::list<int> &P, std::list<std::vector<int> > &pList)
    {
        if(P.size() == (unsigned)Word::get_length())
        {
            std::list<int>::const_iterator Si, Pi;
            std::vector<int> perm(Word::get_length());
            for(Si = S.begin(), Pi = P.begin(); Si != S.end(); Si++, Pi++)
                perm[*Si] = *Pi;
            pList.push_back(perm);
            return;
        }

        Label f, g;
        Signature lf, lg;

        LabelMonomials(ll1, S, f);
        LabelMonomials(ll2, P, g);
        ListSignature(f, S, lf);
        ListSignature(g, P, lg);

        if(EqualSignature(lf, lg))
        {
            std::list<int>::reverse_iterator j;
            Signature::iterator lfi = lf.begin(), lgi = lg.begin();

            if(lf.back().second.size() == 1)
            {
                int num_ins = 0;
                for(; lfi != lf.end(); lfi++, lgi++, num_ins++)
                {
                    S.push_front(lfi->second.front());
                    P.push_front(lgi->second.front());
                }

                FindAllPermutations(ll1, ll2, S, P, pList);
                for(int i = 0; i < num_ins; i++)
                {  S.pop_front();  P.pop_front(); }
            }
            else {

                S.push_front(lfi->second.back());

                for(j = lgi->second.rbegin(); j != lgi->second.rend(); ++j)
                {
                    P.push_front(*j);
                    FindAllPermutations(ll1, ll2, S, P, pList);
                    P.pop_front();
                }

                S.pop_front();
            }
        }
    }

    bool FindPermutation(const Code &c1, const Code &c2, std::list<std::vector<int> > &pList, int level, std::ostream &out)
    {
        std::list<int> S, P;
        std::list<std::set<Word> > ll1, ll2;
        pList.clear();

        c1.lcw_partition_by_levels(ll1, level, hamming_m, degrevlex, false, out);
        c2.lcw_partition_by_levels(ll2, level, hamming_m, degrevlex, false, out);

        clock_t bg_time = clock();
        return FindPermutation(ll1, ll2, S, P, pList, 1, bg_time, out, false);
    }

    void PermutationAnalysis(const Code &cod, const Code &pcod, int level, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out, bool detailed)
    {
        int lv;
        Code perm_cod = cod;
        std::list<std::vector<int> > pList;
        std::list<std::set<Word> > ll1, ll2;
        std::list<std::set<Word> >::iterator i;

        ClearScreen();
        Word::set_wrdout_format(wf);
        out << "\n*** CODE 1 ****\n\n";
        ShowCodeParameters(cod, out);
        cod.lcw_partition_by_levels(ll1, level, wm, mord, false, std::cout);
        out << "\nLeader codewords by levels";
        for(i = ll1.begin(), lv = 1; i != ll1.end(); i++, lv++)
        {
            out << "\nLevel-" << lv << " (" << i->size() << " elements)\n";
            PrintVector(*i, out);
        }

        out << "\n\n\n*** CODE 2 ****\n\n";
        ShowCodeParameters(pcod, out);
        perm_cod.lcw_partition_by_levels(ll2, level, wm, mord, false, std::cout);
        out << "\nLeader codewords by levels";
        for(i = ll2.begin(), lv = 1; i != ll2.end(); i++, lv++)
        {
            out << "\nLevel-" << lv << " (" << i->size() << " elements)\n";
            PrintVector(*i, out);
        }

        out << "\n\nApplying code equivalence algorithm "
            << "to find the permutation.\n"
            << "The set of leader codewords up to level 2 is used as invariant."
            << "\nLeader codewords are computed according "
            << "to ZNTS_AAECC Def. 5 p. 6 using the ";
        ShowGBLAComputingParameters(wm, mord, out);
        out << "\n\n";

        std::list<int> S, P;
        clock_t bg_time = clock();
        if(FindPermutation(ll1, ll2, S, P, pList, 1, bg_time, out, detailed))
            out << "\n*** CODES ARE PERMUTATION EQUIVALENT ***\n";
        else out << "\nCODES ARE NOT PERMUTATION EQUIVALENT !!!\n";

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();
    }

    void AutomorphismGroup(const Code &cod, const Code &pcod, int level, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out)
    {
        int lv, pos;
        Code perm_cod = cod;
        std::list<std::vector<int> > pList;
        std::list<std::set<Word> > ll1, ll2;
        std::list<std::set<Word> >::iterator i;

        ClearScreen();
        Word::set_wrdout_format(wf);
        out << "\n*** CODE 1 ****\n\n";
        ShowCodeParameters(cod, out);
        cod.lcw_partition_by_levels(ll1, level, wm, mord, false, std::cout);
        out << "\nLeader codewords by levels";
        for(i = ll1.begin(), lv = 1; i != ll1.end(); i++, lv++)
        {
            out << "\nLevel-" << lv << " (" << i->size() << " elements)\n";
            PrintVector(*i, out);
        }

        out << "\n\n\n*** CODE 2 ****\n\n";
        ShowCodeParameters(pcod, out);
        perm_cod.lcw_partition_by_levels(ll2, level, wm, mord, false, std::cout);
        out << "\nLeader codewords by levels";
        for(i = ll2.begin(), lv = 1; i != ll2.end(); i++, lv++)
        {
            out << "\nLevel-" << lv << " (" << i->size() << " elements)\n";
            PrintVector(*i, out);
        }

        out << "\n\nApplying code equivalence algorithm "
            << "to find the automorphism group.\n"
            << "The set of leader codewords up to level 2 is used as invariant."
            << "\nLeader codewords are computed according "
            << "to ZNTS_AAECC Def. 5 p. 6 using the ";
        ShowGBLAComputingParameters(wm, mord, out);
        out << "\n\n";

        std::list<int> S, P;
        clock_t bg_time = clock();
        FindAllPermutations(ll1, ll2, S, P, pList);
        if(!pList.empty())
            out << "\n*** CODES ARE PERMUTATION EQUIVALENT ***\n"
                << "Execution time : " << std::setprecision(7)
                << (double)(clock()-bg_time)/CLOCKS_PER_SEC << " seconds.";
        else out << "\nCODES ARE NOT PERMUTATION EQUIVALENT !!!";

        out << "\n\n*** PERMUTATION AUTOMORPHISM GROUP ***\n\n";
        std::list<std::vector<int> >::iterator pi = pList.begin();
        for(pos = 1; pi != pList.end(); pos++, pi++)
        {
            out << pos << ". ";
            PrintPermutationAsCycles(*pi, out);
            if(!IsValidPermutation(ll1, ll2, *pi))
                out << "  <- !! NOT VALID";
            out << "\n";
        }

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();
    }

    void LeaderCodewordsByLevels(const Code &c, word_metric wm,
            monomial_order mord, wrdout_format wf, std::ostream &out, bool std_out)
    {
        std::list<std::set<Word> > ll;
        std::list<std::set<Word> >::iterator i;
        int level = std::numeric_limits<int>::max(), lv;

        ClearScreen();
        if(!std_out)
        {
            ShowCodeParameters(c, std::cout);
            std::cout << "\n\n*** LEADER CODEWORDS BY LEVELS\n"
                      << "\nComputing the leader codewords according "
                      << "to ZNTS_AAECC Def. 5 p. 6 using the ";
            ShowGBLAComputingParameters(wm, mord, std::cout);
            std::cout << "\n\n";
        }

        ShowCodeParameters(c, out);
        out << "\n\n*** LEADER CODEWORDS BY LEVELS\n"
            << "\nComputing the leader codewords according "
            << "to ZNTS_AAECC Def. 5 p. 6 using the ";
        ShowGBLAComputingParameters(wm, mord, out);
        out << "\n\n";

        Word::set_wrdout_format(wf);
        c.lcw_partition_by_levels(ll, level, wm, mord, true, out);

        out << "\n\nLeader codewords by levels";
        for(i = ll.begin(), lv = 1; i != ll.end(); i++, lv++)
        {
            out << "\nLevel-" << lv << " (" << i->size() << " elements)\n";
            PrintVector(*i, out);
        }
        out << "\n\n";

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();

        Word::set_metric(wm);
        Word::set_monomial_order(mord);
    }

    void PermuteCode(Code &c, std::vector<int> &perm,
                     std::ostream &out, bool std_out)
    {
        std::list<std::set<Word> > ll;
        std::list<std::set<Word> >::iterator i;

        ClearScreen();
        ShowCodeParameters(c, out);
        out << "\n\n*** PERMUTED CODE\n\nUsing permutation vector";
        PrintVector(perm, out);
        out << "\nPermutation as cycles: ";
        PrintPermutationAsCycles(perm, out);
        out << "\n\n";

        c.permute(perm);
        ShowCodeParameters(c, out);

        if(!std_out) std::cout << "\n\nSee results in file out.txt in the program directory.\n\n";
        SystemPause();
    }
}
