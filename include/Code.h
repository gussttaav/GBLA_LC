#ifndef CODE_H
#define CODE_H

#include <set>
#include <ctime>
#include <limits>
#include <iomanip>
#include <algorithm>

#include "Word.h"

using namespace GaloisField;

namespace GBLA_LC
{
    struct matphi_elem
    {
        int dom_N, dom_X, img_N;
        bool operator<(const matphi_elem& mp) const
        {
            if(dom_N != mp.dom_N)
                return dom_N < mp.dom_N;
            else return dom_X < mp.dom_X;
        }
    };


    class Code
    {
        public:
            Code();
            Code(const Matrix &checkM);
            Code(const Matrix &genMat, const Matrix &chkMat);
            Code(const Code& cod);

            Code& operator=(const Code& c);

            void random_generate(int len, int dim);
            bool is_codeword(const Word &w) const;
            void permute(const std::vector<int> &perm);
            std::vector<GF> syndrome(const Word &w) const;
            void gen_codewords(std::set<Word> &cwords) const;
            Word encode_vmsg(const std::vector<GF> &mv) const;
            void gen_coset(const Word &w, std::set<Word> &coset) const;
            void gen_codewords_2levels(std::list<std::set<Word> > &cwords) const;
            void generate_variables(std::vector<Word> &X, word_representation wrep);

            void grobner_representation(std::vector<Word> &N, std::set<matphi_elem> &matphi,
                    word_representation wr, word_metric wm, monomial_order mord);
            void leader_codewords(std::set<Word> &lcwords, word_metric wm, monomial_order mord) const;
            void coset_leaders(std::list<std::set<Word> > &leaders, word_metric wm, monomial_order mord);
            void lcw_partition_by_levels(std::list<std::set<Word> > &part,
                    int level, word_metric wm, monomial_order mord, bool print_time, std::ostream &out) const;

            Matrix get_gen_mat() const { return G; }
            Matrix get_check_mat() const { return H;  }
            int get_length() const { return n; }
            int get_dimension() const { return k; }

        private:
            //Code attributes
            int n, k;
            Matrix G, H;

            //Code private methods
            void genmat_from_chkmat();
            void chkmat_from_genmat();
            bool next_mapping(std::vector<int>::iterator fst,
                        std::vector<int>::iterator lst, int fval, int lval) const;

            //GBLA private methods
            Word nextTerm(std::set<Word> &t_list) const;
            int member(const std::vector<GF> &v, const std::list<std::vector<GF> > &G) const;
            bool inCosetLeader(const Word &w, const std::list<std::set<Word> > &leaders) const;
            void insertNexts(const Word &t, std::set<Word> &t_list, const std::vector<Word> &X) const;

            void generate_stdvariables(std::vector<Word> &X);
            void generate_genvariables(std::vector<Word> &X);
            void generate_stdvars_with_exponent(std::vector<Word> &X) const;

            //For debugging purpose
            void print_list(const std::set<Word>& t_list) const;
            void print_irreducibles(int last, const std::vector<Word> &N);
    };

}
#endif // CODE_H
