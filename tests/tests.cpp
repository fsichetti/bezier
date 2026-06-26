#include <catch2/catch_test_macros.hpp>
#include <element_validity.hpp>

using element_validity::fp_t;
using element_validity::Validity;
using element_validity::StaticValidator;
using element_validity::ContinuousValidator;
using element_validity::JacobianEvaluator;

bool almostEq(fp_t a, fp_t b) { return std::fabs(a - b) <= 0.000001; }

TEST_CASE("Standard linear tet validity") {
	StaticValidator<3, 3, 1> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,0.,
 		1.,0.,0.,
 		0.,1.,0.,
 		0.,0.,1.,
	};
	std::vector<int> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::valid);
	JacobianEvaluator<3, 3, 1> evaluator(cp);
	CHECK(evaluator.eval({0,0,0}) == 1);
}

TEST_CASE("Invalid linear tet validity") {
	StaticValidator<3, 3, 1> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		1.,1.,1.,
 		0.,1.,1.,
 		1.,0.,1.,
 		1.,1.,0.,
	};
	std::vector<int> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::invalid);
	JacobianEvaluator<3, 3, 1> evaluator(cp);
	CHECK(evaluator.eval({0,0,0}) == -1.);
}

TEST_CASE("Standard quadratic tet validity") {
	StaticValidator<3, 3, 2> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,0.,
		1.,0.,0.,
		0.,1.,0.,
		0.,0.,1.,
		.5,0.,0.,
		.5,.5,0.,
		0.,.5,0.,
		0.,0.,.5,
		.5,0.,.5,
		0.,.5,.5,
	};
	std::vector<int> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::valid);
	JacobianEvaluator<3, 3, 2> evaluator(cp);
	CHECK(evaluator.eval({0,0,0}) == 1.);
	CHECK(evaluator.eval({1,0,0}) == 1.);
	CHECK(evaluator.eval({0,1,0}) == 1.);
	CHECK(evaluator.eval({0,0,1}) == 1.);
	CHECK(almostEq(evaluator.eval({.3,.3,.3}), 1.));
}


TEST_CASE("Invalid quadratic tet validity") {
	StaticValidator<3, 3, 2> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		.5,.5,.5,
		1.,0.,0.,
		0.,1.,0.,
		0.,0.,1.,
		.5,0.,0.,
		.5,.5,0.,
		0.,.5,0.,
		0.,0.,.5,
		.5,0.,.5,
		0.,.5,.5,
	};
	std::vector<int> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::invalid);
	JacobianEvaluator<3, 3, 2> evaluator(cp);
	CHECK(evaluator.eval({0,0,0}) < 0.);
	CHECK(evaluator.eval({1,0,0}) > .5);
	CHECK(evaluator.eval({0,1,0}) > .5);
	CHECK(evaluator.eval({0,0,1}) > .5);
}

TEST_CASE("Difficult quadratic tet validity") {
	StaticValidator<3, 3, 2> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		-0.42452634383930121231, 	-0.15199946433181563132, 	-0.52560337937709855626,
		-0.44050913176895295508, 	-0.15470197781589292707, 	-0.53435160512323620629,
		-0.42918452277309032671, 	-0.16301150961501897906, 	-0.52689501546960226807,
		-0.43797994733333217576, 	-0.14995527783228912444, 	-0.53289783632873310992,
		-0.43084313006168462045, 	-0.15188149791866326077, 	-0.53189337930183577807,
		-0.43243098252764933242, 	-0.15859376743924027963, 	-0.5309488096272174662,
		-0.42647517075677771059, 	-0.15685809076376799886, 	-0.5264561675879664282,
		-0.4314661156158717814, 	-0.15098007388019624164, 	-0.52468553005828000301,
		-0.43654101656446037127, 	-0.15442916576536502848, 	-0.53114214085428379519,
		-0.43290118494677543026, 	-0.1583252066831398741, 	-0.52536247464088814407,
	};
	std::vector<int> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::uncertain);
}

TEST_CASE("Standard linear tet MTS") {
	ContinuousValidator<3, 3, 1> checker;
	const double prec = .1;
    checker.setPrecisionTarget(prec);
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,	0.,0.,	0.,0.,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
	};
	std::vector<int> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts == 1);
}

TEST_CASE("Standard quadratic tet MTS") {
	ContinuousValidator<3, 3, 2> checker;
	const double prec = .1;
    checker.setPrecisionTarget(prec);
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,	0.,0.,	0.,0.,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
		.5,.5,	0.,0.,	0.,0.,
		.5,.5,	.5,.5,	0.,0.,
		0.,0.,	.5,.5,	0.,0.,
		0.,0.,	0.,0.,	.5,.5,
		.5,.5,	0.,0.,	.5,.5,
		0.,0.,	.5,.5,	.5,.5,
	};
	std::vector<int> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts == 1);
}

TEST_CASE("Invalid quadratic tet MTS") {
	ContinuousValidator<3, 3, 2> checker;
	const double prec = .1;
    checker.setPrecisionTarget(prec);
    const std::vector<double> cp = {
		0.,.5,	0.,.5,	0.,.5,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
		.5,.5,	0.,0.,	0.,0.,
		.5,.5,	.5,.5,	0.,0.,
		0.,0.,	.5,.5,	0.,0.,
		0.,0.,	0.,0.,	.5,.5,
		.5,.5,	0.,0.,	.5,.5,
		0.,0.,	.5,.5,	.5,.5,
	};
	std::vector<int> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts < 1);
}

TEST_CASE("Exact zero at corner quadratic tet MTS") {
	ContinuousValidator<3, 3, 2> checker;
	const double prec = .001;
    checker.setPrecisionTarget(prec);
    const std::vector<double> cp = {
		0.,0.,	0.,0.,	0.,0.,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
		.5,.5,	0.,0.,	0.,0.,
		.5,.25,	.5,.5,	0.,0.,
		0.,0.,	.5,.5,	0.,0.,
		0.,0.,	0.,0.,	.5,.5,
		.5,.5,	0.,0.,	.5,.5,
		0.,0.,	.5,.5,	.5,.5,
	};
	std::vector<int> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts < 1);
	CHECK(mts > 1-prec);
}

TEST_CASE("Standard cubic tetrahedron") {
	StaticValidator<3, 3, 3> checker;
	const std::vector<double> cp = {
		0.0, 0.0, 0.0,
		3.0, 0.0, 0.0,
		0.0, 3.0, 0.0,
		0.0, 0.0, 3.0,
		1.0, 0.0, 0.0,
		2.0, 0.0, 0.0,
		2.0, 1.0, 0.0,
		1.0, 2.0, 0.0,
		0.0, 2.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 0.0, 2.0,
		2.0, 0.0, 1.0,
		1.0, 0.0, 2.0,
		0.0, 2.0, 1.0,
		0.0, 1.0, 2.0,
		1.0, 1.0, 0.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 1.0,
		0.0, 1.0, 1.0,
	};
	std::vector<int> ah;
	const Validity val = checker.isValid(cp, &ah);
	CHECK(val == Validity::valid);
	JacobianEvaluator<3, 3, 3> evaluator(cp);
	CHECK(almostEq(evaluator.eval({0,0,0}), 27.));
	CHECK(almostEq(evaluator.eval({1,0,0}), 27.));
	CHECK(almostEq(evaluator.eval({.1,.1,.1}), 27.));
}
#ifdef INSTANTIATE_3_3_4
TEST_CASE("Standard quartic tetrahedron") {
	StaticValidator<3, 3, 4> checker;
	const std::vector<double> cp = {
		0.0, 0.0, 0.0, 
		4.0, 0.0, 0.0, 
		0.0, 4.0, 0.0, 
		0.0, 0.0, 4.0, 
		1.0, 0.0, 0.0, 
		2.0, 0.0, 0.0, 
		3.0, 0.0, 0.0, 
		3.0, 1.0, 0.0, 
		2.0, 2.0, 0.0, 
		1.0, 3.0, 0.0, 
		0.0, 3.0, 0.0, 
		0.0, 2.0, 0.0, 
		0.0, 1.0, 0.0, 
		0.0, 0.0, 1.0, 
		0.0, 0.0, 2.0, 
		0.0, 0.0, 3.0, 
		3.0, 0.0, 1.0, 
		2.0, 0.0, 2.0, 
		1.0, 0.0, 3.0, 
		0.0, 3.0, 1.0, 
		0.0, 2.0, 2.0, 
		0.0, 1.0, 3.0, 
		1.0, 1.0, 0.0, 
		1.0, 2.0, 0.0, 
		2.0, 1.0, 0.0, 
		1.0, 0.0, 1.0, 
		1.0, 0.0, 2.0, 
		2.0, 0.0, 1.0, 
		2.0, 1.0, 1.0, 
		1.0, 1.0, 2.0, 
		1.0, 2.0, 1.0, 
		0.0, 2.0, 1.0, 
		0.0, 1.0, 2.0, 
		0.0, 1.0, 1.0, 
		1.0, 1.0, 1.0, 
	};
	std::vector<int> ah;
	const Validity val = checker.isValid(cp, &ah);
	CHECK(val == Validity::valid);
	JacobianEvaluator<3, 3, 4> evaluator(cp);
	CHECK(almostEq(evaluator.eval({0,0,0}), 64.));
	CHECK(almostEq(evaluator.eval({1,0,0}), 64.));
	CHECK(almostEq(evaluator.eval({.1,.1,.1}), 64.));
}
#endif
TEST_CASE("Standard linear hex") {
	StaticValidator<3, 1, 1> checker;
	const std::vector<double> cp = {
		0., 0., 0.,
		1., 0., 0.,
		1., 1., 0.,
		0., 1., 0.,
		0., 0., 1.,
		1., 0., 1.,
		1., 1., 1.,
		0., 1., 1.,
	};
	std::vector<int> ah;
	const Validity val = checker.isValid(cp, &ah);
	CHECK(val == Validity::valid);
	JacobianEvaluator<3, 1, 1> evaluator(cp);
	CHECK(almostEq(evaluator.eval({0,0,0}), 1.));
	CHECK(almostEq(evaluator.eval({1,0,0}), 1.));
	CHECK(almostEq(evaluator.eval({.1,.1,.1}), 1.));
}

TEST_CASE("Dynamic standard linear hex") {
	ContinuousValidator<3, 1, 1> checker;
	const std::vector<double> cp = {
		0.,0., 0.,0., 0.,0.,
		1.,1., 0.,0., 0.,0.,
		1.,1., 1.,1., 0.,0.,
		0.,0., 1.,1., 0.,0.,
		0.,0., 0.,0., 1.,1.,
		1.,1., 0.,0., 1.,1.,
		1.,1., 1.,1., 1.,1.,
		0.,0., 1.,1., 1.,1.,
	};
	std::vector<int> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	CHECK(mts == 1.);
}