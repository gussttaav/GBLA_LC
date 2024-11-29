#include "GFElem.h"

namespace GaloisField
{
    //default extension values
    int p = 2,
        q = 2,
        m = 1;

	GFElem alpha(1);

	GFElem::GFElem() : a(m, 0) { }

	GFElem::GFElem(const std::string &str, int posi) : a(m, 0)
	{   parse_from_string(str, posi);	}

	GFElem::GFElem(int elem) : a(m, 0)
	{   a[m-1] = elem%p;	}

	GFElem::GFElem(int *elem) : a(m)
	{
		for (int i = 0; i < m; i++)
			a[i] = elem[i]%p;
	}

	GFElem::GFElem(const GFElem& gfe)
	{   a = gfe.a;	}

	GFElem& GFElem::operator=(const GFElem& gfe)
	{
		if (this == &gfe)
			return *this; // handle self assignment

		a = gfe.a;
		return *this;
	}

	GFElem& GFElem::operator=(const int& gfi)
	{
		for(int i = 0; i < m-1; i++)
			a[i] = 0;
		a[m-1] = gfi%p;
		return *this;
	}

	GFElem GFElem::operator-() const
	{
		GFElem result = *this;

		for (int i = 0; i < m; i++)
			result.a[i] = -a[i];

		return result;
	}

	GFElem& GFElem::operator++()
	{
		a[m-1] = (a[m-1]+1)%p;
		return *this;
	}

	GFElem& GFElem::operator--()
	{
		a[m-1] = (a[m-1]-1)%p;
		return *this;
	}

	int& GFElem::operator[](int index)
	{   return a[index];    }

    const int& GFElem::operator[] (int index) const
    {   return a[index];    }

	GFElem GFElem::operator+(const GFElem& gfe) const
	{
	    if(m == 1) return GFElem(a[0]+gfe.a[0]);

		GFElem result;
		for(int i = 0; i < m; i++)
			result.a[i] = (a[i]+gfe.a[i])%p;
		return result;
	}

	GFElem GFElem::operator+(const int& gfi) const
	{
        GFElem result = *this;
        result.a[m-1] = (result.a[m-1]+gfi)%p;
        return result;
	}

	GFElem& GFElem::operator+=(const GFElem& gfe)
	{   return *this = (*this) + gfe;   }

	GFElem& GFElem::operator+=(const int& gfi)
	{   return *this = (*this) + gfi;     }

	GFElem GFElem::operator-(const GFElem& gfe) const
	{
	    if(m == 1) return GFElem((a[0]-gfe.a[0])%p);

		GFElem result;
		for(int i = 0; i < m; i++)
			result.a[i] = (a[i]-gfe.a[i])%p;
		return result;
	}

	GFElem GFElem::operator-(const int& gfi) const
	{
        GFElem result = *this;
        result.a[m-1] = (result.a[m-1]-gfi)%p;
        return result;
	}

	GFElem& GFElem::operator-=(const GFElem& gfe)
	{   return *this = (*this) - gfe;   }

	GFElem& GFElem::operator-=(const int& gfi)
	{   return *this = (*this) - gfi;     }

	GFElem GFElem::operator*(const GFElem& gfe) const
	{
	    if(m == 1) return GFElem((a[0]*gfe.a[0])%p);

	    //convolution product
	    //prod - product
	    int plen = 2*m-1;
	    int prod[plen], op1[plen], op2[plen], sum;
	    for(int i = 0; i < plen; i++)
        {   prod[i] = 0;  op1[i] = 0;  op2[i] = 0;    }

        for(int i1 = 0, i2 = plen-1; i1 < m; i1++, i2--)
        {
            op1[i1] = a[i1];
            op2[i2] = gfe.a[i1];
        }

	    for(int i = 0, k = plen-1; i < plen; i++, k--)
        {
            sum = 0;
            for(int j = 0; j < plen; j++)
                sum += op1[j]*op2[j];
            //shift left
            for(int j = 0; j < plen-1; j++)
                op2[j] = op2[j+1];
            op2[plen-1] = 0;
            prod[k] = sum%p;
        }

        //Create primitive polynomial from alpha
        int prim_pol[plen], tprod[plen], dprim = m-2, dprod, lt, ltd;
        for(int i = 0; i < plen; i++) prim_pol[i] = 0;
        for(int i = dprim+1, j = 0; i < plen; i++, j++)
            prim_pol[i] = -alpha.a[j];
        prim_pol[dprim] = 1;
        for(dprod = 0; dprod < plen && prod[dprod] == 0; dprod++);

        //division algorithm. prod = prod mod prim_pol
        while(dprod < plen && dprod <= dprim)
        {
            //(lt(prod)/lt(prim_pol))*prim_pol
            lt = prod[dprod];
            ltd = dprim-dprod;
            for(int i = 0; i < dprim; i++) tprod[i] = 0;
            for(int i = dprim; i < plen; i++) tprod[i-ltd] = prim_pol[i];
            for(int i = plen-ltd; i < plen; i++) tprod[i] = 0;

            //prod = prod - (lt(prod)/lt(prim_pol))*prim_pol
            for(int i = dprod; i < plen; i++)
                prod[i] = (prod[i] - lt*tprod[i])%p;
            while(dprod < plen && prod[dprod] == 0) dprod++;
        }

        //copy result
        GFElem result;
        for(int i = plen-1, j = m-1; i >= dprod; i--, j--)
            result.a[j] = prod[i];

        return result;
	}

	GFElem GFElem::operator*(const int& gfi) const
	{
        GFElem result = *this;
        for(int i = 0; i < m; i++)
            result.a[i] = (result.a[i]*gfi)%p;
        return result;
	}

	GFElem& GFElem::operator*=(const GFElem& gfe)
	{   return *this = (*this) * gfe;   }

	GFElem& GFElem::operator*=(const int& gfi)
	{   return *this = (*this) * gfi;     }

	bool GFElem::operator==(const GFElem& gfe) const
	{
	    for(int i = 0; i < m; i++)
            if(a[i] != gfe.a[i])
                return false;
        return true;
	}

	bool GFElem::operator!=(const GFElem& gfe) const
	{   return !((*this) == gfe);   }

	bool GFElem::operator==(const int& gfi) const
	{   return (*this) == GFElem(gfi);  }

    bool GFElem::operator!=(const int& gfi) const
    {   return (*this) != GFElem(gfi);  }

	std::ostream& GFElem::output(std::ostream& out) const
	{
		int i, ua;
		for (i = 0; i < m && a[i] == 0; i++);

		if (i >= m-1)
			out << a[m-1];
		else { // is not a number

			if(a[i] == 1)
				out << "a";
			else if(a[i] == -1)
				out << "-a";
			else out << a[i] << "*a";

			if(i < m-2)
				out << "^" << m-i-1;

			for (int j = i+1; j < m; j++)
			{
				if(a[j] != 0)
				{
					if(a[j] > 0)
						out << " + ";
					else out << " - ";
					ua = a[j] > 0 ? a[j] : -a[j];

					if(ua != 1)
					{
						out << ua;
						if(j < m-2)
							out << "*a^" << m-j-1;
						else if(j == m-2)
							out << "*a";
					}
					else {
						if(j < m-2)
							out << "a^" << m-j-1;
						else if(j == m-2)
							out << "a";
						else out << "1";
					}
				} // end if(a[j] != 0)
			} // end for
		} // end else

		return out;
	}

	std::ostream &operator<<(std::ostream &out, const GFElem& gfe)
    {   return gfe.output(out); }

    bool GFElem::parse_from_string(std::string in_str, int ipos)
    {
        char c = ' ';
        char nbuffer[256];
        std::string snum;
        int ta[m], pos = ipos, st_num, exp,
            sg = 1, val, state = 1, len = in_str.length();
        for(int i = 0; i < m; i++) ta[i] = 0;

        if(ipos >= len) return false;

        while(true)
        {
            if(pos < len)
               c = in_str[pos];

            switch(state)
            {
                case 1: if(pos == len)
                           return false;

                        switch(c)
                        {
                            case '+': sg = 1;
                                      state = 2;
                                break;

                            case '-': sg = -1;
                                      state = 2;
                                break;

                            default: if(isdigit(c))
                                     {
                                         state = 3;
                                         st_num = pos;
                                     }
                                     else if(isalpha(c))
                                     {
                                         val = 1;
                                         state = 4;
                                     }
                                     else if(c != ' ')
                                         return false;
                                break;
                        }
                    break;  //end of case 1

                case 2: if(pos == len)
                           return false;

                        if(isdigit(c))
                        {
                            state = 3;
                            st_num = pos;
                        }
                        else if(isalpha(c))
                        {
                            val = sg;
                            state = 4;
                        }
                        else if(c != ' ')
                            return false;
                    break;

                case 3: while(pos < len && isdigit(in_str[pos])) pos++;
                        snum = in_str.substr(st_num, pos);
                        strcpy(nbuffer, snum.c_str());
                        val = sg*atoi(nbuffer);

                        if(pos == len)
                        {
                            ta[m-1] = (ta[m-1]+val)%p;
                            goto succ_parse;
                        }
                        else{
                            while(pos < len && in_str[pos] == ' ')
                                pos++;
                            if(pos == len)
                            {
                                ta[m-1] = (ta[m-1]+val)%p;
                                goto succ_parse;
                            }
                            c = in_str[pos];

                            switch(c)
                            {
                                case '+': sg = 1;
                                          state = 2;
                                          ta[m-1] = (ta[m-1]+val)%p;
                                    break;

                                case '-': sg = -1;
                                          state = 2;
                                          ta[m-1] = (ta[m-1]+val)%p;
                                    break;

                                case '*': state = 6;
                                    break;

                                default: return false;
                            }
                        }
                    break;

                case 4: if(m == 1) return false;
                        while(pos < len && isalpha(in_str[pos])) pos++;

                        if(pos == len)
                        {
                            ta[m-2] = (ta[m-2]+val)%p;
                            goto succ_parse;
                        }
                        else{
                            while(pos < len && in_str[pos] == ' ')
                                pos++;
                            if(pos == len)
                            {
                                ta[m-2] = (ta[m-2]+val)%p;
                                goto succ_parse;
                            }
                            c = in_str[pos];

                            switch(c)
                            {
                                case '+': sg = 1;
                                          state = 2;
                                          ta[m-2] = (ta[m-2]+val)%p;
                                    break;

                                case '-': sg = -1;
                                          state = 2;
                                          ta[m-2] = (ta[m-2]+val)%p;
                                    break;

                                case '^': state = 8;
                                    break;
                            }
                        }
                    break;

                case 6: if(pos == len)
                           return false;

                        if(isalpha(c))
                            state = 4;
                        else if(c != ' ')
                            return false;
                    break;

                case 8: if(pos == len)
                           return false;

                        if(isdigit(c))
                        {
                            state = 9;
                            st_num = pos;
                        }
                        else if(c != ' ')
                            return false;
                    break;

                case 9: while(pos < len && isdigit(in_str[pos])) pos++;
                        snum = in_str.substr(st_num, pos);
                        strcpy(nbuffer, snum.c_str());
                        exp = atoi(nbuffer);
                        if(exp >= m) return false;

                        if(pos == len)
                        {
                            ta[m-exp-1] = (ta[m-exp-1]+val)%p;
                            goto succ_parse;
                        }
                        else{
                            while(pos < len && in_str[pos] == ' ')
                                pos++;
                            if(pos == len)
                            {
                                ta[m-exp-1] = (ta[m-exp-1]+val)%p;
                                goto succ_parse;
                            }
                            c = in_str[pos];

                            switch(c)
                            {
                                case '+': sg = 1;
                                          state = 2;
                                          ta[m-exp-1] = (ta[m-exp-1]+val)%p;
                                    break;

                                case '-': sg = -1;
                                          state = 2;
                                          ta[m-exp-1] = (ta[m-exp-1]+val)%p;
                                    break;

                                default: return false;
                            }
                        }
                    break;
            }

            pos++;
        }

        //copy elements
        succ_parse:
        for(int i = 0; i < m; i++)
            a[i] = ta[i];

        return true;
    }
}
