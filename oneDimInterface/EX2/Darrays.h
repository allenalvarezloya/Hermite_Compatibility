// Darray classes for arrays of degree 1,2,3, and 4

#ifndef __DARRAYS_H_INCLUDED__

#define __DARRAYS_H_INCLUDED__

class Darray1
{
public:
    Darray1( int nc, int ibeg, int iend);
    Darray1( int ibeg, int iend);
    Darray1();
    ~Darray1() {if( m_data != 0 ) delete[] m_data;}
    void define( int nc, int ibeg, int iend);
    void define( int ibeg, int iend);
    int npts() const {return m_ni;}
    inline double* c_ptr() {return m_data;}
    inline bool in_range( int c, int i)
    {return 1 <= c && c <= m_nc && m_ib <= i && i <= m_ie;}
    // overload parenthesis operator
    inline double& operator()( int c, int i){
        // Turn on array bounds check
#ifdef BZ_DEBUG
        try
        {
            if (!in_range(c,i)) throw 10;
        }
        catch(int e)
        {
            std::cout <<
            "Error Index (c,i) = (" << c << "," <<
            i << ") not in range 1<= c <= " <<
            m_nc << " " << m_ib << " <= i <= " << m_ie<< "\n";
        }
#endif
        return m_data[m_base+m_offc*c+m_offi*i];
    }
    inline double& operator()( int i)
    {
#ifdef BZ_DEBUG
        try{
            if (!in_range(1,i)) throw 10;
        }
        catch(int e)
        {
            std::cout << "Error Index (c,i,j,k) = (" << 1
            << "," << i << ") not in range 1<= c <= "
            << m_nc << " "
            << m_ib << " <= i <= " << m_ie << "\n";
        }
#endif
        return m_data[m_base+m_offi*i+m_offc];
    }
    inline bool is_defined(){return m_data != NULL;}
    static bool m_corder;
    int m_ib, m_ie;
    ssize_t m_base;
    size_t m_offi, m_offc, m_npts;
    void define_offsets();
    void set_value( double scalar );
    void copy( const Darray1& u );
    void writeToFile(char* fileName, int is, int ie);
    int m_nc, m_ni;
private:
    double* m_data;
};

class Darray2
{
public:
    Darray2( int nc, int ibeg, int iend, int jbeg, int jend);
    Darray2( int ibeg, int iend, int jbeg, int jend);
    Darray2();
    ~Darray2() {if( m_data != 0 ) delete[] m_data;}
    void define( int nc, int ibeg, int iend, int jbeg, int jend);
    void define( int ibeg, int iend, int jbeg, int jend);
    inline double* c_ptr() {return m_data;}
    int npts() const {return m_ni*m_nj;}
    inline bool in_range( int c, int i, int j)
    {return 1 <= c && c <= m_nc && m_ib <= i && i <= m_ie && m_jb <= j && j <= m_je;}
    // overload parenthesis operator
    inline double& operator()( int c, int i, int j){
        // Turn on array bounds check
#ifdef BZ_DEBUG
        try
        {
            if (!in_range(c,i,j)) throw 10;
        }
        catch(int e)
        {
            std::cout <<
            "Error Index (c,i,j) = (" << c << "," <<
            i << "," << j <<") not in range 1<= c <= " <<
            m_nc << " " << m_ib << " <= i <= " << m_ie <<
            " " << m_jb << " <= j <= " << m_je<< "\n";
        }
#endif
        return m_data[m_base+m_offc*c+m_offi*i+m_offj*j];
    }
    inline double& operator()( int i, int j)
    {
#ifdef BZ_DEBUG
        try{
            if (!in_range(1,i,j)) throw 10;
        }
        catch(int e)
        {
            std::cout <<
            "Error Index (1,i,j) = (" << 1 << "," <<
            i << "," << j <<") not in range 1<= c <= " <<
            m_nc << " " << m_ib << " <= i <= " << m_ie <<
            " " << m_jb << " <= j <= " << m_je << "\n";
        }
#endif
        return m_data[m_base+m_offi*i+m_offc+m_offj*j];
    }
    inline bool is_defined(){return m_data != NULL;}
    static bool m_corder;
    int m_ib, m_ie;
    int m_jb, m_je;
    ssize_t m_base;
    size_t m_offi, m_offj, m_offc, m_npts;
    void define_offsets();
    void set_value( double scalar );
    void copy( const Darray2& u );
    void writeToFile(char* fileName, int is, int ie, int js, int je);
    int m_nc, m_ni, m_nj;
private:
    double* m_data;
};

class Darray3
{
public:
    Darray3( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend);
    Darray3( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend);
    Darray3();
    ~Darray3() {if( m_data != 0 ) delete[] m_data;}
    void define( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend);
    void define( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend);
    inline double* c_ptr() {return m_data;}
    int npts() const {return m_ni*m_nj*m_nk;}
    inline bool in_range( int c, int i, int j, int k)
    {return 1<= c && c <= m_nc && m_ib <= i && i <= m_ie && m_jb <= j && j <= m_je && m_kb <= k && k <= m_ke;}
    // overload parenthisis operator
    inline double& operator()( int c, int i, int j, int k){
        // Turn on array bounds check
#ifdef BZ_DEBUG
        try
        {
            if (!in_range(c,i,j,k)) throw 10;
        }
        catch(int e)
        {
            std::cout <<
            "Error Index (c,i,j,k) out of bounds" << "\n";
        }
#endif
        return m_data[m_base+m_offc*c+m_offi*i+m_offj*j+m_offk*k];
    }
    inline double& operator()( int i, int j, int k)
    {
#ifdef BZ_DEBUG
        try{
            if (!in_range(1,i,j,k)) throw 10;
        }
        catch(int e)
        {
            std::cout <<
            "Error Index (1,i,j,k) out of bounds" << "\n";
        }
#endif
        return m_data[m_base+m_offi*i+m_offc+m_offj*j+m_offk*k];
    }
    inline bool is_defined(){return m_data != NULL;}
    static bool m_corder;
    int m_ib, m_ie;
    int m_jb, m_je;
    int m_kb, m_ke;
    ssize_t m_base;
    size_t m_offi, m_offj, m_offk, m_offc, m_npts;
    void define_offsets();
    void set_value( double scalar );
    void copy( const Darray3& u );
    int m_nc, m_ni, m_nj, m_nk;
private:
    double* m_data;
    
};




class Darray4
{
public:
    Darray4( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int lbeg, int lend);
    Darray4( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int lbeg, int lend);
    Darray4();
    ~Darray4() {if( m_data != 0 ) delete[] m_data;}
    void define( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int lbeg, int lend);
    void define( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int lbeg, int lend);
    inline double* c_ptr() {return m_data;}
    int npts() const {return m_ni*m_nj*m_nk*m_nl;}
    inline bool in_range( int c, int i, int j, int k, int l)
    {return 1<= c && c <= m_nc && m_ib <= i && i <= m_ie && m_jb <= j && j <= m_je && m_kb <= k && k <= m_ke && m_lb <= l && l <= m_le;}
    // overload parenthisis operator
    inline double& operator()( int c, int i, int j, int k, int l){
        // Turn on array bounds check
#ifdef BZ_DEBUG
        try
        {
            if (!in_range(c,i,j,k,l)) throw 10;
        }
        catch(int e)
        {
            std::cout <<
            "Error Index (c,i,j,k,l) out of bounds" << "\n";
        }
#endif
        return m_data[m_base+m_offc*c+m_offi*i+m_offj*j+m_offk*k+m_offl*l];
    }
    inline double& operator()( int i, int j, int k, int l)
    {
#ifdef BZ_DEBUG
        try{
            if (!in_range(1,i,j,k,l)) throw 10;
        }
        catch(int e)
        {
            std::cout <<
            "Error Index (1,i,j,k) out of bounds" << "\n";
        }
#endif
        return m_data[m_base+m_offi*i+m_offc+m_offj*j+m_offk*k+m_offl*l];
    }
    inline bool is_defined(){return m_data != NULL;}
    static bool m_corder;
    int m_ib, m_ie;
    int m_jb, m_je;
    int m_kb, m_ke;
    int m_lb, m_le;
    ssize_t m_base;
    size_t m_offi, m_offj, m_offk, m_offl,  m_offc, m_npts;
    void define_offsets();
    void set_value( double scalar );
    void copy( const Darray4& u );
    int m_nc, m_ni, m_nj, m_nk, m_nl;
    double getMax();
private:
    double* m_data;
    
};

#endif


