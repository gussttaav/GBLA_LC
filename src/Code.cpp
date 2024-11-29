
#include "Code.h"

namespace GBLA_LC
{

    Code::Code(): n(0), k(0), G(), H() { }

    Code::Code(const Matrix &checkM) : n(checkM.get_cols()),
                    k(n-checkM.get_rows()), G(k, n), H(checkM)
    {
        genmat_from_chkmat();
        Word::set_length(n);
    }

    Code::Code(const Matrix &genMat, const Matrix &chkMat) : n(genMat.get_cols()),
                    k(genMat.get_rows()), G(genMat), H(chkMat)
    {
        Word::set_length(n);
    }



    Code::Code(const Code& cod)
    {
        n = cod.n;
        k = cod.k;
        H = cod.H;
        G = cod.G;
        Word::set_length(n);
    }

    Code& Code::operator=(const Code& c)
    {
        if(this == &c)
            return *this;

        n = c.n;
        k = c.k;
        H = c.H;
        G = c.G;
        Word::set_length(n);

        return *this;
    }

    void Code::random_generate(int len, int dim)
    {
        assert(len > 0 && dim > 0 && len > dim);

        //(I, A) I-(k*k) Identity matrix

        n = len;
        k = dim;
        G = Matrix(k, n);
        Word::set_length(n);

        for(int i = 0; i < k; i++)
            G(i,i) = 1;

        srand(time(NULL));
        for(int i = 0; i < k; i++)
            for(int j = k; j < n; j++)
                G(i,j) = GF(rand() % q);
        chkmat_from_genmat();

        //generate random permutation
        int rvalue, pos = 0;
        std::vector<int> perm(n, -1);
        while(pos < n)
        {
            rvalue = rand() % n;
            if(std::find(perm.begin(), perm.end(), rvalue) == perm.end())
            { perm[pos] = rvalue; pos++; }
        }

        permute(perm);
    }

    void Code::grobner_representation(std::vector<Word> &N, std::set<matphi_elem> &matphi,
                    word_representation wr, word_metric wm, monomial_order mord)
    {
        if(H.is_empty()) return;

        Word w;
        std::vector<GF> v;
        std::set<Word> List;
        std::vector<Word> X;
        std::list<std::vector<GF> > G;
        std::vector<Word>::size_type g;
        int j, r = -1;

        Word::set_metric(wm);
        Word::set_monomial_order(mord);
        generate_variables(X, wr);
        List.insert(Word());
        N.resize(pow(q, n-k));
        matphi_elem mphie;

        //int y = 1;
        //Word::set_wrdout_format(wf_monomial);
        while(!List.empty())
        {
            //std::cout << "Iteracion " << y << "\nList:";
            //print_list(List);

            w = nextTerm(List);
            v = syndrome(w);
            j = member(v, G);

            //std::cout << "w = " << w << "\n\n";
            //y++;
            if(j != -1)
            {
                mphi_upd:
                for(int i = 0; i <= r; i++)
                    for(g = 0; g < X.size(); g++)
                    {
                        if(w == N[i]*X[g])
                        {
                            mphie.dom_N = i;
                            mphie.dom_X = g;
                            mphie.img_N = j;
                            matphi.insert(mphie);
                        }
                    }
            }
            else{
                r = r+1;
                N[r] = w;
                G.push_back(v);
                insertNexts(w, List, X);

                //Update matphi
                j = r;
                goto mphi_upd;
            }
        }
    }

    void Code::coset_leaders(std::list<std::set<Word> > &leaders,
                             word_metric wm, monomial_order mord = degrevlex)
    {
        if(H.is_empty()) return;

        Word w;
        int j, r = -1;
        std::vector<GF> v;
        std::vector<Word> X, N;
        std::set<Word> List, nset;
        std::list<std::vector<GF> > G;
        std::list<std::set<Word> >::iterator itj;
        std::list<std::vector<GF> >::iterator itg;
        std::pair<std::set<Word>::iterator, bool> ret;

        Word::set_metric(wm);
        Word::set_monomial_order(mord);
        leaders.clear();
        List.insert(Word());
        N.resize(pow(q, n-k));
        generate_stdvars_with_exponent(X); //SF

        while(!List.empty())
        {
            w = nextTerm(List);
            v = syndrome(w);
            j = member(v, G);

            if(j != -1)
            {
                if(w.get_weight() <= N[j].get_weight())
                {
                    itj = leaders.begin();
                    std::advance(itj, j);

                    if(w.get_weight() < N[j].get_weight())
                    {
                        N[j] = w;
                        itg = G.begin();
                        std::advance(itg, j);
                        (*itg) = v;
                        (*itj).clear();
                    }

                    ret = (*itj).insert(w);
                    if(ret.second)
                        insertNexts(w, List, X);
                }

                //Inserting in List according to Criterion 2 in ZNTS_AAECC Def. 4, p. 5.
                //**Assume standard representation
                if(w.get_weight() == N[j].get_weight()+1 && j > 0)
                {
                    for(int i = 0; i < n; i++)
                    {
                        if(w[i] != 0)
                        {
                            Word cl = w;
                            cl[i] = 0;
                            cl.apply_inv_morphism();
                            cl.calculate_weigth();

                            if(!inCosetLeader(cl, leaders))
                            {
                                v = syndrome(cl);
                                j = member(v, G);

                                if((j != -1) && (cl.get_weight() == N[j].get_weight()))
                                {
                                    itj = leaders.begin();
                                    std::advance(itj, j);
                                    (*itj).insert(cl);
                                    insertNexts(cl, List, X);
                                }
                            }
                        }
                    }
                } //end criterion 2
            }
            else{
                r = r+1;
                N[r] = w;
                nset.clear();
                G.push_back(v);
                nset.insert(w);
                leaders.push_back(nset);
                insertNexts(w, List, X);
            }
        }
    }

    void Code::leader_codewords(std::set<Word> &lcwords, word_metric wm, monomial_order mord) const
    {
        if(H.is_empty()) return;

        Word w, t;
        int j, r = -1;
        std::vector<GF> v;
        std::vector<Word> X, N;
        std::set<Word> List, nset;
        std::list<std::vector<GF> > G;
        std::list<std::set<Word> > CL;
        std::list<std::set<Word> >::iterator itj;
        std::list<std::vector<GF> >::iterator itg;
        std::pair<std::set<Word>::iterator, bool> ret;

        Word::set_metric(wm);
        Word::set_monomial_order(mord);
        CL.clear();
        lcwords.clear();
        List.insert(Word());
        N.resize(pow(q, n-k));
        generate_stdvars_with_exponent(X); //SF

        while(!List.empty())
        {
            w = nextTerm(List);
            v = syndrome(w);
            j = member(v, G);

            if(j != -1)
            {
                itj = CL.begin();
                std::advance(itj, j);

                if(w.get_weight() <= N[j].get_weight())
                {
                    if(w.get_weight() < N[j].get_weight())
                    {
                        N[j] = w;
                        itg = G.begin();
                        std::advance(itg, j);
                        (*itg) = v;
                        (*itj).clear();
                    }

                    ret = (*itj).insert(w);
                    if(ret.second)
                        insertNexts(w, List, X);
                }

                //Inserting in List according to Criterion 2 in ZNTS_AAECC Def. 4, p. 5.
                //**Assume standard representation
                if(w.get_weight() == N[j].get_weight()+1 && j > 0)
                {
                    for(int i = 0; i < n; i++)
                    {
                        if(w[i] != 0)
                        {
                            Word cl = w;
                            cl[i] = 0;

                            if(!inCosetLeader(cl, CL))
                            {
                                v = syndrome(cl);
                                j = member(v, G);

                                cl.calculate_weigth();
                                if((j != -1) && (cl.get_weight() == N[j].get_weight()))
                                {
                                    itj = CL.begin();
                                    std::advance(itj, j);
                                    cl.apply_inv_morphism();
                                    ret = (*itj).insert(cl);
                                    if(ret.second)
                                        insertNexts(cl, List, X);
                                }
                            }
                        }
                    }
                } //end criterion 2

                //Leader codewords computation according to ZNTS_AAECC Def. 5 p. 6
                //**Assume standard representation
                for(std::set<Word>::iterator s = (*itj).begin(); s != (*itj).end(); s++)
                    if(w != (*s))
                    {
                        for(int i = 0; i < n; i++)
                            t[i] = w[i]-(*s)[i];

                        t.apply_inv_morphism();
                        t.calculate_weigth();
                        ret = lcwords.insert(t);

                        bool match_found = false;
                        if(ret.second) //leader codeword added
                        {
                            GF gf_tmp;
                            GFElem fp;
                            int ei, ej, tmp;
                            //remove if w = n1 + eij and n1 is not a coset leader
                            for(ei = 0; ei < n && !match_found; ei++)
                            {
                                if(w[ei] != 0)
                                {
                                    fp = w[ei].as_polynomial();
                                    for(ej = 0; ej < q && !match_found; ej++)
                                        if(fp[ej] != 0)
                                        {
                                            tmp = fp[ej];  gf_tmp = w[ei];
                                            fp[ej] = 0;    w[ei] = fp;
                                            match_found = inCosetLeader(w, CL);
                                            fp[ej] = tmp;    w[ei] = gf_tmp;
                                        }
                                }
                            }
                            if(!match_found)
                                lcwords.erase(ret.first);
                        }
                    }
            }
            else{
                r = r+1;
                N[r] = w;
                nset.clear();
                G.push_back(v);
                nset.insert(w);
                CL.push_back(nset);
                insertNexts(w, List, X);
            }
        }
    }


    void Code::lcw_partition_by_levels(std::list<std::set<Word> > &part,
            int level, word_metric wm, monomial_order mord, bool print_time, std::ostream &out) const
    {
        if(H.is_empty()) return;

        GF gf_tmp;
        GFElem fp;
        Word w, t;
        std::vector<GF> v;
        std::vector<Word> X, N;
        std::set<Word> List, nset;
        std::list<std::vector<GF> > G;
        std::list<std::set<Word> > CL;
        std::list<std::vector<GF> >::iterator itg;
        std::list<std::set<Word> >::iterator itj, li;
        std::pair<std::set<Word>::iterator, bool> ret;
        int j, r = -1, aw, wt, lw = 0, ei, ej, tmp/*, lc_num = 0, step = 15*/;
        bool match_found, part_found;

        Word::set_metric(wm);
        Word::set_monomial_order(mord);
        CL.clear();
        part.clear();
        List.insert(Word());
        N.resize(pow(q, n-k));
        generate_stdvars_with_exponent(X); //SF

        clock_t bg_time = clock();
        //lw = std::numeric_limits<int>::max();

        do{
            w = nextTerm(List);
            v = syndrome(w);
            j = member(v, G);

            aw = w.get_weight();

            if(j != -1)
            {
                itj = CL.begin();
                std::advance(itj, j);

                if(w.get_weight() <= N[j].get_weight())
                {
                    if(w.get_weight() < N[j].get_weight())
                    {
                        N[j] = w;
                        itg = G.begin();
                        std::advance(itg, j);
                        (*itg) = v;
                        (*itj).clear();
                    }

                    ret = (*itj).insert(w);
                    if(ret.second)
                        insertNexts(w, List, X);
                }

                //Inserting in List according to Criterion 2 in ZNTS_AAECC Def. 4, p. 5.
                //**Assume standard representation
                if(w.get_weight() == N[j].get_weight()+1 && j > 0)
                {
                    for(int i = 0; i < n; i++)
                    {
                        if(w[i] != 0)
                        {
                            Word cl = w;
                            cl[i] = 0;

                            if(!inCosetLeader(cl, CL))
                            {
                                v = syndrome(cl);
                                j = member(v, G);

                                cl.calculate_weigth();
                                if((j != -1) && (cl.get_weight() == N[j].get_weight()))
                                {
                                    itj = CL.begin();
                                    std::advance(itj, j);
                                    cl.apply_inv_morphism();
                                    ret = (*itj).insert(cl);
                                    if(ret.second)
                                        insertNexts(cl, List, X);
                                }
                            }
                        }
                    }
                } //end criterion 2

                //Leader codewords computation according to ZNTS_AAECC Def. 5 p. 6
                //**Assume standard representation
                for(std::set<Word>::iterator s = (*itj).begin(); s != (*itj).end(); s++)
                    if(w != (*s))
                    {
                        for(int i = 0; i < n; i++)
                            t[i] = w[i]-(*s)[i];

                        t.apply_inv_morphism();
                        t.calculate_weigth();
                        wt = t.get_weight();

                        //create partition by Levels with leader codeword t
                        part_found = false;
                        for(li = part.begin(); li != part.end() && wt >= li->begin()->get_weight(); li++)
                            if(li->begin()->get_weight() == wt)
                            {  part_found = true;  break;  }

                        if(part_found)
                            ret = li->insert(t);
                        else{
                            std::set<Word> npart;
                            ret = npart.insert(t);
                            if(part.size() < (unsigned)level)
                            {
                                part.insert(li, npart);
                                if(part.size() == (unsigned)level)
                                    lw = part.back().begin()->get_weight();

                                if(print_time)
                                   out << "Partition of weigth " << t.get_weight() << " created. Elapsed time: "
                                       << std::setprecision(7) << (double)(clock()-bg_time)/CLOCKS_PER_SEC << " seconds.\n";
                            }
                            else if(part.size() == (unsigned)level && li != part.end())
                            {
                                li->clear();
                                ret = li->insert(t);
                                lw = part.back().begin()->get_weight();

                                if(print_time)
                                   out << "Partition of weigth " << t.get_weight() << " created. Elapsed time: "
                                       << std::setprecision(7) << (double)(clock()-bg_time)/CLOCKS_PER_SEC << " seconds.\n";
                            }
                        }

                        if(ret.second) //leader codeword added
                        {
                            //remove if w = n1 + eij and n1 is not a coset leader
                            match_found = false;
                            for(ei = 0; ei < n && !match_found; ei++)
                            {
                                if(w[ei] != 0)
                                {
                                    fp = w[ei].as_polynomial();
                                    for(ej = 0; ej < q && !match_found; ej++)
                                        if(fp[ej] != 0)
                                        {
                                            tmp = fp[ej];  gf_tmp = w[ei];
                                            fp[ej] = 0;    w[ei] = fp;
                                            match_found = inCosetLeader(w, CL);
                                            fp[ej] = tmp;    w[ei] = gf_tmp;
                                        }
                                }
                            }
                            if(!match_found)
                                li->erase(ret.first);
                        }
                    }
            }
            else{
                r = r+1;
                N[r] = w;
                nset.clear();
                G.push_back(v);
                nset.insert(w);
                CL.push_back(nset);
                insertNexts(w, List, X);
            }

        }while(!List.empty() && ((2*aw-1) <= lw || lw == 0));

        if(print_time)
         out << "Overall time: " << std::setprecision(7) << (double)(clock()-bg_time)/CLOCKS_PER_SEC << " seconds.\n";
    }


    bool Code::inCosetLeader(const Word &w, const std::list<std::set<Word> > &leaders) const
    {
        std::set<Word>::const_iterator its;
        std::list<std::set<Word> >::const_iterator itl;

        for(itl = leaders.begin(); itl != leaders.end(); itl++)
        {
            its = std::find((*itl).begin(), (*itl).end(), w);
            if(its != (*itl).end())
                return true;
        }

        return false;
    }

    Word Code::nextTerm(std::set<Word> &t_list) const
    {
        Word fst_elem = *(t_list.begin());
        t_list.erase(t_list.begin());
        return fst_elem;
    }

    std::vector<GF> Code::syndrome(const Word &w) const
    {
        GF sum;
        int slen = n-k;
        std::vector<GF> synd(slen);

        for(std::vector<GF>::size_type i = 0; i < synd.size(); i++)
        {
            sum = 0;
            for(int j = 0; j < n; j++)
                sum += H(i,j) * w[j];
            synd[i] = sum;
        }

        return synd;
    }

    int Code::member(const std::vector<GF> &v, const std::list<std::vector<GF> > &G) const
    {
        std::list<std::vector<GF> >::const_iterator it = std::find(G.begin(), G.end(), v);
        return std::find(G.begin(), G.end(), v) != G.end() ? std::distance(G.begin(), it) : -1;
    }

    void Code::insertNexts(const Word &t, std::set<Word> &t_list, const std::vector<Word> &X) const
    {
        for(std::vector<Word>::size_type i = 0; i < X.size(); i++)
            t_list.insert(X[i]*t);
    }

    //assume check matrix in standard form
    void Code::genmat_from_chkmat()
    {
        if(H.is_empty()) return;

        //Create Identity first
        for(int i = 0; i < k; i++)
            G(i,i) = 1;

        for(int i = 0; i < n-k; i++)
            for(int j = 0; j < k; j++)
               G(j,k+i) = -H(i,j);
    }

    //assume generator matrix in standard form
    void Code::chkmat_from_genmat()
    {
        H = Matrix(n-k, n);

        for(int i = 0, i1 = k; i < n-k; i++, i1++)
            for(int j = 0; j < k; j++)
               H(i,j) = -G(j,i1);

        for(int i = 0, j = k; j < n; i++, j++)
            H(i,j) = 1;
    }

    void Code::generate_variables(std::vector<Word> &X,
                        word_representation wrep = standard_rep)
    {
        Word::set_representation(wrep);

        switch(wrep)
        {
            case standard_rep:
                    generate_stdvariables(X);
                break;

            case general_rep:
                    generate_genvariables(X);
                break;

            default:
                    generate_stdvariables(X);
                break;
        }
    }

    void Code::generate_stdvariables(std::vector<Word> &X)
    {
        X.clear();
        X.resize(n*m);
        Word::set_representation(standard_rep);

        for(std::vector<Word>::size_type i = 0; i < X.size(); i++)
            X[i].set_monomial_var_exponent(i, 1);
    }

    void Code::generate_genvariables(std::vector<Word> &X)
    {
        //not implemented yet
    }

    void Code::generate_stdvars_with_exponent(std::vector<Word> &X) const
    {
        X.clear();
        X.resize(n*(q-1));

        //Generate all variables in standard form
        for(std::vector<Word>::size_type i = 0, h = 0; i < X.size(); h++)
            for(int j = 1; j < q; j++, i++)
            {
                X[i][h] = GF(j);
                X[i].apply_inv_morphism();
            }
    }


    Word Code::encode_vmsg(const std::vector<GF> &mv) const
    {   return Word(mv*G);  }


    bool Code::next_mapping(std::vector<int>::iterator fst,
                    std::vector<int>::iterator lst, int fval, int lval) const
    {
        if(lst == fst)
            return false;

        do{
            if(++(*(--lst)) != lval)
                return true;
            *lst = fval;
        } while(lst != fst);

        return false;
    }


    //only for two levels
    void Code::gen_codewords_2levels(std::list<std::set<Word> > &cwords) const
    {
        Word cw;
        cwords.clear();
        std::set<Word> nset;
        std::vector<GF> v(k);
        std::vector<int> vi(k, 0);
        int wt, lw[2] = {-1, -1};

        //a list with two empty set
        cwords.push_back(nset);
        cwords.push_back(nset);
        std::list<std::set<Word> >::iterator li;
        next_mapping(vi.begin(), vi.end(), 0, q);

        //generate codewords
        do{
            for(int i = 0; i < k; i++)
                v[i] = GF(vi[i]);

            cw = encode_vmsg(v);
            wt = cw.get_weight();
            li = cwords.begin();

            if(lw[0] == -1 || lw[1] == -1)
            {
                if(lw[0] == -1)
                {
                    lw[0] = wt;
                    li->insert(cw);
                }
                else
                {
                    if(lw[0] == wt)
                        li->insert(cw);
                    else{
                        li++;
                        lw[1] = wt;
                        li->insert(cw);
                    }
                }
            }
            else{
                if(wt == lw[0])
                    li->insert(cw);
                else if(wt == lw[1])
                {   li++; li->insert(cw); }
                else{
                    if(wt < lw[0] && wt < lw[1])
                    {
                        if(lw[0] < lw[1])
                        {   li++; li->clear(); li->insert(cw); lw[1] = wt; }
                        else
                        {   li->clear(); li->insert(cw); lw[0] = wt; }
                    }
                    else if(wt < lw[0])
                    {   li->clear(); li->insert(cw); lw[0] = wt; }
                    else if(wt < lw[1])
                    {   li++; li->clear(); li->insert(cw); lw[1] = wt; }
                }
            }

        }
        while(next_mapping(vi.begin(), vi.end(), 0, q));

    }

    void Code::gen_codewords(std::set<Word> &cwords) const
    {
        cwords.clear();
        std::vector<GF> v(k);
        std::vector<int> vi(k, 0);

        do{
            for(int i = 0; i < k; i++)
                v[i] = GF(vi[i]);
            cwords.insert(encode_vmsg(v));
        }
        while(next_mapping(vi.begin(), vi.end(), 0, q));
    }

    void Code::gen_coset(const Word &w, std::set<Word> &coset) const
    {
        coset.clear();
        std::vector<GF> v(k);
        std::vector<int> vi(k, 0);

        do{
            for(int i = 0; i < k; i++)
                v[i] = GF(vi[i]);
            coset.insert(w*encode_vmsg(v));
        }
        while(next_mapping(vi.begin(), vi.end(), 0, q));
    }


    void Code::permute(const std::vector<int> &perm)
    {
        G.permute_cols(perm);
        H.permute_cols(perm);
    }

    bool Code::is_codeword(const Word &w) const
    {
        std::vector<GF> v = syndrome(w);
        std::vector<GF>::iterator i = v.begin();

        for(; i != v.end(); i++)
            if(*i != 0) return false;

        return true;
    }

    //For debugging purpose
    void Code::print_list(const std::set<Word>& t_list) const
    {
        std::set<Word>::const_iterator x;
        for(x = t_list.begin(); x != t_list.end(); x++)
            std::cout << *x << "  ";
        std::cout << std::endl;
    }

    void Code::print_irreducibles(int last, const std::vector<Word> &N)
    {
        for(int i = 0; i <= last; i++)
            std::cout << N[i] << " ";

        std::cout << std::endl;
    }

}
