"""New groups from set operations on existing groups: intersection, union,
difference, complement, and the opposite set of any result."""
import unittest

import numpy as np

from tests.fixtures import make_sample, make_view
from tests.test_clustering_markers import api_client

# A 4 x 4 grid; g0 values are 1, 3, 5, ... 31 (row-major) and g1 values 2, 4, ... 32
A = [dict(column="g0", op=">", value="10", logic="AND")]     # g0 > 10: observations 5..15
B = [dict(column="g0", op="<", value="20", logic="AND")]     # g0 < 20: observations 0..9


class CoreSetOperationTests(unittest.TestCase):
    def setUp(self):
        self.groups = make_view(kind="bins").groups
        for name, conditions in (("A", A), ("B", B)):
            prepared = self.groups.prepare_condition_group(
                name, [(c["column"], c["op"], c["value"], c["logic"]) for c in conditions])
            self.groups.commit_condition_group(*prepared)
        self.a, self.b = self.groups.membership("A"), self.groups.membership("B")

    def test_each_operation_and_its_opposite(self):
        expected = {"intersection": self.a & self.b, "union": self.a | self.b, "difference": self.a & ~self.b}
        symbols = {"intersection": "A ∩ B", "union": "A ∪ B", "difference": "A − B"}
        for operation, mask in expected.items():
            with self.subTest(operation=operation):
                result, description = self.groups.combine(operation, "A", "B")
                np.testing.assert_array_equal(result, mask)
                self.assertEqual(description, symbols[operation])
                opposite, description = self.groups.combine(operation, "A", "B", opposite=True)
                np.testing.assert_array_equal(opposite, ~mask)
                self.assertEqual(description, f"not ({symbols[operation]})")
        complement, description = self.groups.combine("complement", "A")
        np.testing.assert_array_equal(complement, ~self.a)
        self.assertEqual(description, "not A")
        self.assertEqual(int(self.a.sum() + complement.sum()), len(self.a))   # within the observations

    def test_invalid_requests(self):
        for args, message in ((("xor", "A", "B"), "Choose one of"), (("union", "A"), "second group"),
                              (("union", "A", "A"), "two different"), (("union", "A", "nope"), "Unknown group")):
            with self.subTest(args=args), self.assertRaisesRegex(ValueError, message):
                self.groups.combine(*args)


class ApiSetOperationTests(unittest.TestCase):
    def setUp(self):
        self.client, self.session = api_client(self, make_sample("bins"))
        for name, conditions in (("A", A), ("B", B)):
            self.assertEqual(self.client.post("/api/layers/bins/groups",
                                              json=dict(name=name, conditions=conditions)).status_code, 200)

    def combine(self, **body):
        return self.client.post("/api/layers/bins/groups/combine", json=body)

    def test_preview_then_create_groups_with_fixed_membership(self):
        preview = self.combine(operation="intersection", first="A", second="B", preview=True).json()
        self.assertEqual((preview["size"], preview["description"], preview["total"]), (5, "A ∩ B", 16))
        self.assertNotIn("A ∩ B", [g["name"] for g in self.client.get("/api/layers/bins/groups").json()["groups"]])

        created = self.combine(operation="intersection", first="A", second="B").json()
        self.assertEqual((created["group"], created["size"]), ("A ∩ B", 5))
        group = next(g for g in created["groups"] if g["name"] == "A ∩ B")
        self.assertEqual((group["size"], group["description"]), (5, "A ∩ B"))

        named = self.combine(operation="difference", first="A", second="B", name="Only high").json()
        self.assertEqual((named["group"], named["size"]), ("Only high", 6))           # observations 10..15
        self.assertEqual(self.combine(operation="union", first="A", second="B", name="Only high").status_code, 422)
        self.assertEqual(self.combine(operation="complement", first="Only high").json()["size"], 10)
        # The opposite of the union is empty: refused rather than creating an empty group
        empty = self.combine(operation="union", first="A", second="B", opposite=True)
        self.assertEqual(empty.status_code, 422)
        self.assertIn("is empty", empty.json()["detail"])

        # Membership is fixed at creation: deleting an operand keeps the result
        self.client.delete("/api/layers/bins/groups/B")
        groups = {g["name"]: g["size"] for g in self.client.get("/api/layers/bins/groups").json()["groups"]}
        self.assertEqual(groups["A ∩ B"], 5)

    def test_result_groups_are_offered_to_de(self):
        self.combine(operation="union", first="A", second="B")
        options = self.client.get("/api/layers/bins/de/options").json()
        self.assertIn("A ∪ B", str(options))

    def test_unknown_operation_or_group(self):
        self.assertEqual(self.combine(operation="xor", first="A", second="B").status_code, 422)
        self.assertEqual(self.combine(operation="union", first="A", second="nope").status_code, 422)


if __name__ == "__main__":
    unittest.main()
