#define CONFIG 512, 16, uint16_t, 62, uint64_t

#include <iostream>
#include "nfl/prng/FastGaussianNoise.hpp"
#include <nfl/poly_p.hpp>
#include "tools.h"
#define SIGMA 4
#define Berr 1 // From XPIR 5*security_bits/2;
#define A_bits 17	// AbsBitPerCipher / polyDegree
const uint64_t bitmask = (1ULL << A_bits)-1; 	

typedef  nfl::poly_from_modulus<uint64_t, 512, 62> poly;
typedef  nfl::poly_from_modulus<uint16_t, 512, 16> clearpoly;

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
	//   b = (a*s) % f + e * A + m;
	b = nfl::non_uniform(Berr);
	// Adjustments and addition to plaintext
	for(int i=0;i<b.degree;i++) {
		b(0,i) = (( b(0,i) << (A_bits+1) ) + message(0,i));
		if(b(0,i)>b.get_modulus(0)) b(0,i)-=b.get_modulus(0);
	}	
	b.ntt_pow_phi();
	a = nfl::uniform();	
	b = b + nfl::shoup(a * s, sprime);
}

template <class P, class Q>
__attribute__((noinline)) static void encryptNFL(P& a, P& b, Q const & message, P const & s, P const & sprime, nfl::FastGaussianNoise<uint16_t, typename P::value_type, 2> *g_prng)
{
	P tmp;
	plonge(tmp,message);
	encryptNFL(a,b,tmp,s,sprime,g_prng);
}

template <class P>
	__attribute__((noinline)) static void decryptNFL(P& tmp, P const & a, P const& b, P const& s, P const& sprime)
	{
		const uint64_t magicConst = (1ULL<<(s.nbits+1))-a.get_modulus(0);
		tmp = b - nfl::shoup(a * s, sprime);
		tmp.invntt_pow_invphi();
		for(auto & v : tmp)
		{
			//v = (v>a.get_modulus(0)/2) ? ((v + magicConst) & bitmask) : (v & bitmask);
			v = (v % a.get_modulus(0))& bitmask;
		}
	}
	
template <class P, class Q>
	__attribute__((noinline)) static void decryptNFL(Q& res, P const & a, P const& b, P const& s, P const& sprime)
	{
		P tmp;
		decryptNFL(tmp,a,b,s,sprime);
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
		double nbr_bits=floor(( (p_size - 1) - log2(nbsum) - log2(Berr) -log2(nbmul)) / 2.0);
		std::cout << "Berr " << Berr << " nbsums "<<nbsum<< " p_size "<<p_size<<" nmul "<< nbmul << " nbr_bits " << nbr_bits<<std::endl;
		std::cout << "Modulus (ciph/clear) : " << ciph.get_modulus(0)<< " / " << clear.get_modulus(0)<<std::endl;
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
    encryptNFL(ca, cb, p, s, sprime, &g_prng);
	// for (int i=0;i<degree-1;i++) {ea(0,i)=ca(0,i)*ca(0,i);eb(0,i)=cb(0,i)*cb(0,i);}
	ea = ca * ca ;
	eb = cb * cb;
    decryptNFL(result, ea, eb, s, sprime);
	std::cout << "dec(enc(p2*p2)=" << result << std::endl;
	
	
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

