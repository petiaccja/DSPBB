#include <catch2/catch.hpp>
#include <dspbb/Utility/Interval.hpp>

using namespace dspbb;



TEST_CASE("Interval positive offset", "[Interval]") {
	Interval i{ 1, 2 };
	i += 2;
	REQUIRE(i.first == 3);
	REQUIRE(i.last == 4);
}

TEST_CASE("Interval negative offset", "[Interval]") {
	Interval i{ 1, 2 };
	i -= 2;
	REQUIRE(i.first == -1);
	REQUIRE(i.last == 0);
}

TEST_CASE("Interval disjoint", "[Interval]") {
	Interval i1{ 1, 2 };
	Interval i2{ 3,6 };
	Interval i3{ 5,7 };
	REQUIRE(IsDisjoint(i1, i2));
	REQUIRE(!IsDisjoint(i2, i3));
	REQUIRE(IsDisjoint(i1, i3));
}

TEST_CASE("Interval intersection", "[Interval]") {
	Interval i1{ 1, 2 };
	Interval i2{ 3,6 };
	Interval i3{ 5,7 };
	REQUIRE(Intersection(i1, i2).Size() == 0);
	REQUIRE(Intersection(i2, i3) == Interval{5,6});
	REQUIRE(Intersection(i1, i3).Size() == 0);
}

TEST_CASE("Interval encompassing union", "[Interval]") {
	Interval i1{ 1, 2 };
	Interval i2{ 3,6 };
	Interval i3{ 5,7 };
	REQUIRE(EncompassingUnion(i1, i2) == Interval(1,6));
	REQUIRE(EncompassingUnion(i2, i3) == Interval{ 3, 7 });
	REQUIRE(EncompassingUnion(i1, i3) == Interval(1,7));
}

TEST_CASE("Interval union", "[Interval]") {
	Interval i1{ 1, 2 };
	Interval i2{ 3,6 };
	Interval i3{ 5,7 };
	REQUIRE(!Union(i1, i2));
	REQUIRE(Union(i2, i3).value() == Interval{ 3, 7 });
	REQUIRE(!Union(i1, i3));
}
