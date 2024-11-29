
#include "GBLA_LC.h"

#ifdef _WIN32
    #define CLEAR_CMD "cls"
    #define PAUSE_CMD "PAUSE"
#else
    #define CLEAR_CMD "clear"
    #define PAUSE_CMD "read -p 'Press any key to continue...'"
#endif

#define EXIT_PROG   16

using namespace GBLA_LC;

//primality test by brute force
//works fine for small numbers used in this program
bool is_prime(int num)
{
    if(num < 2) return false;

    for(int i = 2; i < num; i++)
    {
        if(num%i == 0)
            return false;
    }
    return true;
}


void CreateCode(Code &cod)
{
    Matrix g, h;
    char dec = '1';
    int ford, clen, cdim;
    std::stringstream ss;
    std::ifstream cfile("codes.txt");
    std::string st, cname, alpha = "1";
    ClearScreen();

    while(dec != '5')
    {
        std::cout << "Create Linear Code\n"
                  << "  1. From file \"codes.txt\"\n"
                  << "  2. Randomly\n"
                  << "  3. Manually\n"
                  << "  4. Back\n\n"
                  << "Select an option number: ";
        std::cin >> dec;
        ClearScreen();

        switch(dec)
        {
            case '1': std::cout << "Loading linear code from file\n"
                                << "Enter the code name: ";
                      std::cin >> cname;
                      while(!ReadCodeFromFile(cfile, cname, cod) && dec != 'n')
                      {
                         std::cout << "ERROR: Linear code " << cname << " could not be "
                                   << "load from file \"codes.txt\". Try again (y/n) ? ";
                         std::cin >> dec;
                      }
                      if(dec == '1') dec = '5';
                break;

            case '2': std::cout << "Generating linear code randomly\n\n"
                                << "Field order: ";
                      std::cin >> ford;
                      if(!is_prime(ford))
                      {
                          std::cout << "Field primitive element: ";
                          std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                          getline(std::cin, alpha);
                      }
                      std::cout << "Code length: ";
                      std::cin >> clen;
                      std::cout << "Code dimension: ";
                      std::cin >> cdim;

                      GF::set_order(ford, alpha);
                      cod.random_generate(clen, cdim);
                      dec = '5';
                break;

            case '3': std::cout << "Creating linear code manually\n\n"
                                << "Field order: ";
                      std::cin >> ford;
                      if(!is_prime(ford))
                      {
                          std::cout << "Field primitive element: ";
                          std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                          getline(std::cin, alpha);
                      }
                      std::cout << "Code length: ";
                      std::cin >> clen;
                      std::cout << "Code dimension: ";
                      std::cin >> cdim;
                      GF::set_order(ford, alpha);

                      std::cout << "Generating Matrix: ";
                      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                      getline(std::cin, st);
                      g = ReadGFMatrixFromStr(st, cdim, clen);
                      std::cout << "Check Matrix: ";
                      getline(std::cin, st);
                      h = ReadGFMatrixFromStr(st, clen-cdim, clen);
                      cod = Code(g, h);
                      dec = '5';
                break;


            case '4': dec = '5';
        }
    }
}

void CreateCodeForEquivalenceTest(Code &cod, Code &pcod)
{
    Matrix g, h;
    char dec = '1';
    int clen = cod.get_length(), cdim = cod.get_dimension();
    std::ifstream cfile("codes.txt");
    std::string st, cname;
    ClearScreen();

    while(dec != '5')
    {
        std::cout << "Create Linear Code for Equivalence Test\n"
                  << "  1. From file \"codes.txt\"\n"
                  << "  2. Randomly\n"
                  << "  3. Manually\n"
                  << "  4. Back\n\n"
                  << "Select an option number: ";
        std::cin >> dec;
        ClearScreen();

        switch(dec)
        {
            case '1': std::cout << "Loading linear code from file\n"
                                << "Enter the code name: ";
                      std::cin >> cname;
                      while(!ReadCodeFromFile(cfile, cname, pcod) && dec != 'n')
                      {
                         std::cout << "ERROR: Linear code " << cname << " could not be "
                                   << "load from file \"codes.txt\". Try again (y/n) ? ";
                         std::cin >> dec;
                      }
                      if(dec == '1') dec = '5';
                break;

            case '2': std::cout << "Generating linear code randomly\n\n";
                      pcod.random_generate(clen, cdim);
                      dec = '5';
                break;

            case '3': std::cout << "Creating linear code manually\n\n";
                      std::cout << "Generating Matrix: ";
                      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                      getline(std::cin, st);
                      g = ReadGFMatrixFromStr(st, cdim, clen);
                      std::cout << "Check Matrix: ";
                      getline(std::cin, st);
                      h = ReadGFMatrixFromStr(st, clen-cdim, clen);
                      pcod = Code(g, h);
                      dec = '5';
                break;


            case '4': dec = '5';
        }
    }
}


word_metric SelectWordMetric()
{
    int op;
    word_metric wm;
    ClearScreen();

    std::cout << "Select word metric\n"
              << "  1. Hamming\n"
              << "  2. Rank\n"
              << "Choose an option number: ";
    std::cin >> op;
    wm = op == 1 ? hamming_m : rank_m;

    return wm;
}

monomial_order SelectMonomialOrder()
{
    int op;
    monomial_order mord;
    ClearScreen();

    std::cout << "Select monomial order\n"
              << "  1. Lexicographical (lex)\n"
              << "  2. Degree lexicographical (deglex)\n"
              << "  3. Degree reverse lexicographical (degrevlex)\n"
              << "Choose an option number: ";
    std::cin >> op;

    switch(op)
    {
        case 1: mord = lex; break;
        case 2: mord = deglex; break;
        case 3: mord = degrevlex; break;
        default: mord = lex;
    }

    return mord;
}

wrdout_format SelectWordFormat()
{
    int op;
    wrdout_format wf;
    ClearScreen();

    std::cout << "Select word output format\n"
              << "  1. Vector\n"
              << "  2. Monomial\n"
              << "  3. Vector:monomial\n"
              << "  4. Vector standard\n"
              << "Choose an option number: ";
    std::cin >> op;

    switch(op)
    {
        case 1: wf = wf_vector; break;
        case 2: wf = wf_monomial; break;
        case 3: wf = wf_pair; break;
        case 4: wf = wf_monomial_sv; break;
        default: wf = wf_vector;
    }

    return wf;
}

bool SelectOutput()
{
    int op;
    ClearScreen();

    std::cout << "Select where to show the program results\n"
              << "  1. Screen\n"
              << "  2. File (out.txt by default in program root directory)\n"
              << "Choose an option number: ";
    std::cin >> op;

    return op == 1 ? true : false;
}


void ReadPermutationVector(const Code &c, std::vector<int> &perm)
{
    char dec;
    std::string st;
    int pos = 0, value, psize = c.get_length();
    perm.clear(); perm.resize(psize, -1);

    ClearScreen();
    std::cout << "[" << c.get_length() << ", " << c.get_dimension()
              << "] Linear code over GF(" << GaloisField::q << ")\n"
              << "Generator Matrix: \n" << c.get_gen_mat()
              << "\nCheck Matrix: \n" << c.get_check_mat() << "\n\n"
              << "Permutation vector randomly (y/n): ";
    std::cin >> dec;
    if(dec != 'y'){
        std::cout << "\n\nType an integer permutation vector of length " << psize << "\n"
                  << "Integers separated by space. Ex. 1 2 3 -> is a vector of length 3.\n\n";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        getline(std::cin, st);
        std::stringstream ss(st);
        int k = 0;
        while(k < psize && ss >> value)
        { perm[k] = value; k++; }
        std::cout << "\n\nInserted vector\n";
    }
    else{
        srand(time(NULL));
        while(pos < psize)
        {
            value = rand() % psize;
            if(std::find(perm.begin(), perm.end(), value) == perm.end())
            { perm[pos] = value; pos++; }
        }

        std::cout << "\n\nRandom permutation generated\n";
    }

    std::vector<int>::iterator i = perm.begin();
    std::cout << "[" << *i;
    for(i++; i != perm.end(); i++)
        std::cout << ", " << *i;
    std::cout << "]\n\n";

    SystemPause();
}

int SelectFromMenu(const Code &c, word_metric wm,
        monomial_order mo, wrdout_format wf, bool std_out)
{
    int dec;
    ClearScreen();
    std::cout << "Computational Algebra Group (2018)\n"
              << "Mathematics Department - University of Oriente - Cuba\n\n";
    if(c.get_length() != 0)
    {
        std::cout << "[" << c.get_length() << ", " << c.get_dimension()
                  << "] Linear code over GF(" << GaloisField::q << ")\n"
                  << "Generator Matrix: \n" << c.get_gen_mat()
                  << "\nCheck Matrix: \n" << c.get_check_mat() << "\n\n";
    }
    std::cout << "MAIN MENU\n"
              << "--------------------------------------------\n"
              << "  1. Generate a new linear code\n"
              << "  2. Select word metric ";
    if(wm == hamming_m)
        std::cout << "(hamming)\n";
    else std::cout << "(rank)\n";
    std::cout << "  3. Select monomial order ";
    switch(mo){
        case lex: std::cout << "(lex)\n"; break;
        case deglex: std::cout << "(deglex)\n"; break;
        case degrevlex: std::cout << "(degrevlex)\n"; break;
    }
    std::cout << "  4. Select word output format ";
    switch(wf){
        case wf_vector: std::cout << "(vector)\n"; break;
        case wf_monomial: std::cout << "(monomial)\n"; break;
        case wf_pair: std::cout << "(vector:monomial)\n"; break;
        case wf_monomial_sv: std::cout << "(vector std)\n"; break;
    }
    std::cout << "--------------------------------------------\n"
              << "Warning! - A code must be generated first for using these functions\n"
              << "  5. Grobner Representation\n"
              << "  6. Coset leaders\n"
              << "  7. Coset leaders analysis\n"
              << "  8. Codewords\n"
              << "  9. Leader codewords\n"
              << "  10. Leader codewords by levels\n"
              << "  11. Permute\n"
              << "  12. Permutation equivalence test\n"
              << "  13. Detailed permutation equivalence algorithm\n"
              << "  14. Automorphism group\n"
              << "--------------------------------------------\n"
              << "  15. Results output ";
    if(std_out)
        std::cout << "(screen)\n";
    else std::cout << "(file)\n";
    std::cout << "  16. Exit\n\n"
              << "Select an option number: ";
    std::cin >> dec;

    return dec;
}


int main()
{
    std::ofstream file;
    std::ostream *outf = &std::cout;

    int op;
    Code cod, pcod;
    bool std_out = true;
    std::vector<int> perm;
    wrdout_format wf = wf_monomial;
    word_metric wm = hamming_m;
    monomial_order mord = degrevlex;

    while((op=SelectFromMenu(cod, wm, mord, wf, std_out)) != EXIT_PROG)
    {
        if(!std_out && !file.is_open())
            file.open("out.txt", std::ios::out);

        switch(op)
        {
            case 1: CreateCode(cod); break;
            case 2: wm = SelectWordMetric(); break;
            case 3: mord = SelectMonomialOrder(); break;
            case 4: wf = SelectWordFormat(); break;
            case 5: GrobnerRepresentation(cod, standard_rep, wm, mord, wf, *outf, std_out); break;
            case 6: CosetLeaders(cod, wm, mord, wf, *outf, std_out); break;
            case 7: TestCosetLeaders(cod, wm, mord, wf, *outf, std_out); break;
            case 8: Codewords(cod, wf, *outf, std_out); break;
            case 9: LeaderCodewords(cod, wm, mord, wf, *outf, std_out); break;
            case 10: LeaderCodewordsByLevels(cod, wm, mord, wf, *outf, std_out); break;
            case 11: ReadPermutationVector(cod, perm);
                     PermuteCode(cod, perm, *outf, std_out);
                break;
            case 12: CreateCodeForEquivalenceTest(cod, pcod);
                     PermutationAnalysis(cod, pcod, 2, wm, mord, wf, *outf, std_out, false);
                break;
            case 13: CreateCodeForEquivalenceTest(cod, pcod);
                     PermutationAnalysis(cod, pcod, 2, wm, mord, wf, *outf, std_out, true);
                break;
            case 14: CreateCodeForEquivalenceTest(cod, pcod);
                     AutomorphismGroup(cod, pcod, 2, wm, mord, wf, *outf, std_out);
                break;
            case 15: std_out = SelectOutput();
                     if(std_out) outf = &std::cout;
                     else outf = &file;
                break;
        }

        if(file.is_open())
            file.close();
    }

    //--GF::set_order(3);
    //Matrix G = ReadGFMatrixFromStr(genMat, 7, 15);
    //Matrix H = ReadGFMatrixFromStr(chkMat, 8, 15);

    //Code cod(G, H);

    //--Code cod;
    //--cod.random_generate(32, 20);

    /*for(int i = 0, len = 11, dim = 8; len <= 15; i++, len++, dim++)
        for(int j = 0; j < 5; j++, i++)
        {
            std::cout << "======== RANDOM EXPERIMENT #" << i+1 << " ========\n\n";
            cod.random_generate(len, dim);
            PermutationAnalysis(cod, 2);
            std::cout << "\n\n\n";
        }*/

    //ReadCodeFromFile(cfile, code_num, cod);
    //--PermutationAnalysis(cod, 2);
    //LeaderCodewordsAnalysisByLevels(cod, std::numeric_limits<int>::max(), hamming_m, degrevlex, wf_vector);

    //FindPermutation(cod, perm_cod, perm);
    //**GrobnerRepresentation(cod, standard_rep, hamming_m, degrevlex, wf_monomial);
    //TestCosetLeaders(cod, rank_m, degrevlex, wf_monomial_sv);
    //CosetLeaders(cod, rank_m, degrevlex, wf_monomial_sv);
    //CosetLeaders(cod, hamming_m, degrevlex, wf_vector);
    //CompareCosetLeadersGAP(cod, "C2-10-4cls.txt", hamming_m, degrevlex, wf_monomial_sv);
    //LeaderCodewords(cod, hamming_m, degrevlex, wf_monomial_sv);

    GF::cleanup_optables();
    //cfile.close();

    /*GF::set_order(9, "a+1");*/

    /*std::cout << "Field elements\n";
    GF::show_field_elems();
    std::cout << "\nADD TABLE\n";
    GF::show_add_table();
    std::cout << "\nMULTIPLICATION TABLE\n";
    GF::show_mult_table();
    std::cout << "\n";*/


   return 0;
}
