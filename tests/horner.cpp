#define CONFIG 1024, 16, uint16_t, 62, uint64_t

#include <iostream>
#include "nfl/prng/FastGaussianNoise.hpp"
#include <nfl/poly_p.hpp>
#include "tools.h"
#include <bitset>

#define SIGMA 4
#define Berr 200 // From XPIR 5*security_bits/2;
#define A_bits 16	// AbsBitPerCipher / polyDegree
const uint64_t bitmask = (1ULL << A_bits)-1; 	
const uint64_t antibitmask = 0x7FFFFFFFFFFFFFFF - bitmask ; 	

typedef  nfl::poly_from_modulus<uint64_t, 1024, 62> poly;
typedef  nfl::poly_from_modulus<uint16_t, 1024, 16> clearpoly;

template <class P, class Q>
static void plonge(P& out, Q& in) {
	for(int j=0;j<out.nmoduli;j++) {
		for(int i=0;i<out.degree;i++) {
			out(j,i)=in(j,i);
		}
	}
}

template <class P>
__attribute__((noinline)) static void encryptNFL(P& a, P& b, P const & message, P const & s, P const & sprime, nfl::FastGaussianNoise<uint16_t, typename P::value_type, 2> *g_prng)
{
	//   b = a.s + e.Amp + m;
	a=nfl::uniform();
	P e = nfl::non_uniform(2*Berr-1);
	for(auto & v : e) {
		v = ((v  << A_bits) );
		if(v>a.get_modulus(0)) v=a.get_modulus(0)-v;
	}
	e = e + message;
	e.ntt_pow_phi();
	b = nfl::shoup(a * s, sprime) + e;
}

template <class P, class Q>
__attribute__((noinline)) static void encryptNFL(P& a, P& b, Q const & message, P const & s, P const & sprime, nfl::FastGaussianNoise<uint16_t, typename P::value_type, 2> *g_prng)
{
	P tmp;
	plonge(tmp,message);
	//tmp.ntt_pow_phi();
	encryptNFL(a,b,tmp,s,sprime,g_prng);
}

template <class P>
	__attribute__((noinline)) static void decryptNFL(P& tmp, P const & a, P const& b, P const& s, P const& sprime)
	{
		// m = b - a.s -e.Amp
		// const uint64_t magicConst = (1ULL<<(s.nbits+1))-a.get_modulus(0);
		tmp = (b - nfl::shoup(a * s, sprime));// % a.get_modulus(0);
		tmp.invntt_pow_invphi();
		for(auto & v : tmp)
		{
			//v = (v>a.get_modulus(0)/2) ? ((v + magicConst) & bitmask) : (v & bitmask);
//			v = v % a.get_modulus(0);
			v = v & bitmask;
		}
	}
	
template <class P, class Q>
	__attribute__((noinline)) static void decryptNFL(Q& res, P const & a, P const& b, P const& s, P const& sprime)
	{
		P tmp;
		decryptNFL(tmp,a,b,s,sprime);
		//tmp.invntt_pow_invphi();
		plonge(res,tmp);
	}
	
template <class P>
	__attribute__((noinline)) static void keyGenNFL(P& s, P& sprime, nfl::FastGaussianNoise<uint16_t, typename P::value_type, 2> *g_prng )
	{
		s = *alloc_aligned<poly, 32>(1, nfl::gaussian<uint16_t, typename P::value_type, 2>(g_prng));
		sprime = *alloc_aligned<poly, 32>(1);
		s.ntt_pow_phi();
		sprime = nfl::compute_shoup(s);
	}
	
template <class P, class Q>
	__attribute__((noinline)) static void computeNoise(P& ciph, Q& clear)
	{
		double nbsum=ciph.degree;
		double p_size=ciph.aggregated_modulus_bit_size;
		double nbmul=ciph.degree;
		double nbr_bits=floor(( (p_size - 1) - log2(nbsum) - log2(2*Berr-1) -log2(nbmul)) / 2.0);
		std::cout << "Berr " << Berr << " nbsums "<<nbsum<< " p_size "<<p_size<<" nmul "<< nbmul << " nbr_bits " << nbr_bits<<std::endl;
		std::cout << "Modulus (ciph/clear) : " << ciph.get_modulus(0)<< " / " << clear.get_modulus(0)<<std::endl;
		std::cout << "Bitmasks (plain/anti) : " << bitmask<< " / " << antibitmask <<std::endl;
	}

template <size_t degree, size_t modulus_clear, class T_clear, size_t modulus_cipher, class T_cipher>
bool run() {

	clearpoly p{1,2,3};
	T_clear x0=3, px0;
	
	// On cherche à évaluer p en x0 soit p(x0)

	/******************* En clair *****************************/
	
	for(int i=degree-1;i>=0;i=i-1) {
		px0=px0*x0+p(0,i);
	}
	assert(px0==34);
	
	/******************* En chiffré maintenant *****************************/
	poly s, sprime;
	clearpoly result;
	poly ca, cb, da, db, ea, eb, res;
	
	// This step generates a secret key
	nfl::FastGaussianNoise<T_clear, T_cipher, 2> g_prng(SIGMA, 128, 1<<10);
	keyGenNFL(s,sprime,&g_prng);
	computeNoise(s, p);
	
	// Test enc/dec result=Dec(Enc(p));
    encryptNFL(ca, cb, p, s, sprime, &g_prng);
    decryptNFL(result, ca, cb, s, sprime);
	bool testEncDec = true; for (int i=0;i<degree-1;i++) testEncDec &= (result(0,i)==p(0,i));
	std::cout << "Test dec(enc(p))==p : " << ((testEncDec) ? " OK" : " KO") << std::endl;
	
	// Test HFE add result=Dec(Enc(p)+Enc(p))
    encryptNFL(ca, cb, p, s, sprime, &g_prng);
	ea = ca + ca ; eb = cb + cb;
    decryptNFL(result, ea, eb, s, sprime);
	bool testFHEAdd = true; for (int i=0;i<degree-1;i++) testFHEAdd &= (result(0,i)==p(0,i)*2);
	std::cout << "Test dec(enc(p)enc(p)) : " << ((testFHEAdd) ? " OK" : " KO") << std::endl;
	
	// Test HFE mul
	p={1,2,3};//p.ntt_pow_phi();
    encryptNFL(ca, cb, p, s, sprime, &g_prng);
	 for (int i=0;i<degree-1;i++) {ea(0,i)=ca(0,i)*ca(0,i);eb(0,i)=cb(0,i)*cb(0,i);}
	// ea = ca * ca ;
	// eb = cb * cb;
    decryptNFL(result, ea, eb, s, sprime);
	std::cout << "dec(enc(p2*p2)=" << result << std::endl;
	clearpoly real_p2{1,4,10,12,9};
	bool testFHEMul = true; for (int i=0;i<degree-1;i++) testFHEMul &= (result(0,i)==real_p2(0,i));
	std::cout << "Test dec(enc(p))==p : " << ((testFHEMul) ? " OK" : " KO") << std::endl;
	
	//
	// uint64_t px0=3,accumul=1;
	// poly ra;
	// poly rb;
	// for(int i=0;i<degree;i=i+1) {
	//  		ra(0,i)=accumul*ca(0,i);
	//  		rb(0,i)=accumul*cb(0,i);
	// 	accumul=accumul*px0;
	//  	}
	// ra = ra*px0+ca;
// 	rb = rb*px0+cb;

	//  decryptNFL(res, ra, rb, s, sprime);
	// uint64_t resultat=0;
	// for(auto & v : res) resultat += v;
	//
	// std::cout << "dec(enc(p2)=" << res << std::endl;
	// std::cout << "p(x0)="<<(resultat )<< std::endl;
//
// 	//T_clear x0=3;
//
// 	// On cherche à évaluer p en x0 soit p(x0)
//
// 	acc=0;
// 	for(int i=degree-1;i>=0;i=i-1) {
// 		acc=acc*x0+p(0,i);
// 	}
//
// 	assert(acc==34);
	
	return true;
}

int main(int argc, char const* argv[]) {
	return not run<CONFIG>();
}

