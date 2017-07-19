#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions.hpp>
#include <complex>
#include <iostream>

/*
	Only change PRECISION
*/

const long long PRECISION = 100, INACC = 4;

typedef boost::multiprecision::number<
	boost::multiprecision::cpp_dec_float<PRECISION + INACC> > arb;

arb alladi_grinstead()
{
	arb c;

	for(unsigned long long n = 1;; ++n) {
		const arb last = c;
		c += (boost::math::zeta<arb>(n + 1) - 1) / n;
		if(c == last) break;
	}

	return exp(c - 1);
}

arb aperys()
{
	return boost::math::zeta<arb>(3);
}

arb buffons()
{
	return 2 / boost::math::constants::pi<arb>();
}

arb catalans()
{
	const arb PI = boost::math::constants::pi<arb>();
	return 0.125 * (boost::math::trigamma<arb>(0.25) - PI * PI);
}

arb champernowne()
{
	arb c;
	for(unsigned long long n = 1;; ++n) {
		arb sub;

		for(unsigned long long k = 1; k <= n; ++k) {
			sub += floor(log10(k));
		}

		const arb last = c;
		c += n / pow(10, n + sub);
		if(c == last) break;
	}

	return c;
}

arb delian()
{
	return boost::math::cbrt<arb>(2);
}

arb dottie()
{
	arb x;
	std::string precomp, postcomp;

	for(x = 1; x== 1 || precomp != postcomp;) {
		precomp = static_cast<std::string>(x);
		precomp.resize(PRECISION);

		x -= (cos(x) - x) / (-sin(x) - 1);

		postcomp = static_cast<std::string>(x);
		postcomp.resize(PRECISION);
	}

	return x;
}

arb e()
{
	return boost::math::constants::e<arb>();
}

arb erdos_borwein()
{
	arb e;

	for(unsigned long long n = 1;; ++n) {
		const arb last = e;
		e += 1 / (pow(static_cast<arb>(2), n) - 1);
		if(e == last) break;
	}

	return e;
}

arb euler_mascheroni()
{
	return boost::math::constants::euler<arb>();
}

arb favard(const unsigned long long r = 0)
{
	if(r % 2 == 0) {
		return (-4 / boost::math::constants::pi<arb>()) *
			(pow(-2, static_cast<arb>(-2) * (r + 1)) / boost::math::tgamma<arb>(r + 1)) *
			(boost::math::polygamma<arb>(r, static_cast<arb>(0.25)) - boost::math::polygamma<arb>(r, static_cast<arb>(0.75)));
	} else {
		return (4 / boost::math::constants::pi<arb>()) *
			((1 - pow(2, -(static_cast<arb>(r) + 1))) * boost::math::zeta<arb>(r + 1));
	}
}

arb gauss()
{
	const arb ROOT_TWO = boost::math::constants::root_two<arb>(),
		PI = boost::math::constants::pi<arb>();
	return pow(boost::math::tgamma<arb>(0.25), 2)
		/ (2 * ROOT_TWO * pow(PI, 3 / static_cast<arb>(2)));
}

arb gelfond_schneider()
{
	return pow(2, boost::math::constants::root_two<arb>());
}

arb gelfonds()
{
	return pow(boost::math::constants::e<arb>(), boost::math::constants::pi<arb>());
}

arb giesekings()
{
	return (9 - boost::math::trigamma<arb>(2 / static_cast<arb>(3)) +
		boost::math::trigamma<arb>(4 / static_cast<arb>(3))) /
		(4 * boost::math::constants::root_three<arb>());
}

arb golden_ratio()
{
	return boost::math::constants::phi<arb>();
}

std::complex<arb> i()
{
	return std::complex<arb>(0,1);
}

arb inverse_golden_ratio()
{
	return boost::math::constants::phi<arb>() - 1;
}

arb khinchin_levy()
{
	return pow(boost::math::constants::pi<arb>(), 2) / (12 * log(static_cast<arb>(2)));
}

arb knuth()
{
	return (1 - (1 / boost::math::constants::root_three<arb>())) / 2;
}

arb levys()
{
	return exp(pow(boost::math::constants::pi<arb>(), 2) / (12 * log(static_cast<arb>(2))));
}

arb liebs()
{
	return (8 * boost::math::constants::root_three<arb>()) / 9;
}

arb lochs()
{
	return (6 * log(static_cast<arb>(2)) * log(static_cast<arb>(10))) /
		pow(boost::math::constants::pi<arb>(), 2);
}

arb nortons()
{
	const arb PI = boost::math::constants::pi<arb>(),
		EULER = boost::math::constants::euler<arb>(),
		GLAISHER = boost::math::constants::glaisher<arb>(),
		PI_SQR = pow(PI, 2),
		LOG_TWO = log(static_cast<arb>(2));
	return -((PI_SQR - 6 * LOG_TWO * (-3 + 2 * EULER + LOG_TWO + 24 * log(GLAISHER) - 2 * log(PI))) / PI_SQR);
}

arb omega()
{
	arb omega;
	std::string precomp, postcomp;

	for(omega = 0; omega == 0 || precomp != postcomp;) {
		precomp = static_cast<std::string>(omega);
		precomp.resize(PRECISION);

		omega -= ((omega * exp(omega)) - 1) /
			(exp(omega) * (omega + 1) - ((omega + 2) * (omega * exp(omega) - 1) / ((2 * omega) + 2)));

		postcomp = static_cast<std::string>(omega);
		postcomp.resize(PRECISION);
	}

	return omega;
}

arb one()
{
	return 1;
}

arb pi()
{
	return boost::math::constants::pi<arb>();
}

arb plastic_number()
{
	return (boost::math::cbrt<arb>(108 + 12 * sqrt(static_cast<arb>(69))) + 
		boost::math::cbrt<arb>(108 - 12 * sqrt(static_cast<arb>(69)))) / static_cast<arb>(6);
}

arb pogsons()
{
	return pow(100, 1 / static_cast<arb>(5));
}

arb polyas_random_walk()
{
	const arb PI = boost::math::constants::pi<arb>(),
		PI_CBD = pow(PI, 3),
		ROOT_SIX = sqrt(static_cast<arb>(6));
	return 1 - 1/((ROOT_SIX / (32 * PI_CBD)) *
        boost::math::tgamma<arb>(1 / static_cast<arb>(24)) *
        boost::math::tgamma<arb>(5 / static_cast<arb>(24)) *
        boost::math::tgamma<arb>(7 / static_cast<arb>(24)) *
        boost::math::tgamma<arb>(11 / static_cast<arb>(24)));
}

arb porters()
{
	const arb PI = boost::math::constants::pi<arb>(),
		GLAISHER = boost::math::constants::glaisher<arb>();
	return ((6 * log(static_cast<arb>(2)) * (48 * log(GLAISHER) - log(static_cast<arb>(2)) - 4 * log(PI) - 2))
		/ pow(PI, 2)) - (1 / static_cast<arb>(2));
}

arb prince_ruperts_cube()
{
	return (3 * boost::math::constants::root_two<arb>()) / 4;
}

arb pythagoras()
{
	return boost::math::constants::root_two<arb>();
}

arb robbins()
{
	const arb PI = boost::math::constants::pi<arb>(),
		ROOT_TWO = boost::math::constants::root_two<arb>(),
		ROOT_THREE = boost::math::constants::root_three<arb>();
	return ((4 + 17 * ROOT_TWO - 6 * ROOT_THREE - 7 * PI) / 105)
		+ (log(1 + ROOT_TWO) / 5) + ((2 * log(2 + ROOT_THREE)) / 5);
}

arb sierpinski_k()
{
	const arb PI = boost::math::constants::pi<arb>(),
		E = boost::math::constants::e<arb>(),
		EULER = boost::math::constants::euler<arb>();
	return PI * log((4 * pow(PI, 3) * pow(E, 2 * EULER)) / pow(boost::math::tgamma<arb>(0.25), 4));
}

arb sierpinski_s()
{
	const arb PI = boost::math::constants::pi<arb>(),
		E = boost::math::constants::e<arb>(),
		EULER = boost::math::constants::euler<arb>();
	return log((4 * pow(PI, 3) * pow(E, 2 * EULER)) / pow(boost::math::tgamma<arb>(0.25), 4));
}

arb silver_ratio()
{
	return boost::math::constants::root_two<arb>() + 1;
}

arb theodorus()
{
	return boost::math::constants::root_three<arb>();
}

arb twenty_vertex_entropy()
{
	return (3 * boost::math::constants::root_three<arb>()) / 2;
}

arb weierstrass()
{
	const arb PI = boost::math::constants::pi<arb>(),
		E = boost::math::constants::e<arb>();
	return (pow(2, static_cast<arb>(1.25)) * sqrt(PI) * pow(E, PI / 8)) /
		pow(boost::math::tgamma<arb>(0.25), 2);
}

arb wylers() {
	const arb PI = boost::math::constants::pi<arb>();
	return (9 / (8 * pow(PI, 4))) * pow(pow(PI, 5) / 1920, 0.25);
}

arb zero()
{
	return 0;
}

int main()
{
	std::cout << std::fixed << std::setprecision(PRECISION)
		<< favard() << '\n';
	return 0;
}
