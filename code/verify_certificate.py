"""Independent checker of the emitted order-six algebraic certificate.

This file imports only the Python standard library, never the generator or its
helpers. It reads the proposed polynomial identities and verifies every defining
equation, every recursive ordered-tree substitution, every exact-flow value,
and the small rank certificates. It does not solve for polynomial formulas.

Run: python verify_certificate.py
"""


def require(condition, detail="verification obligation failed"):
    """A proof obligation that remains active under python -O."""
    if not condition:
        raise ValueError(detail)


import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from math import factorial
from pathlib import Path

Q = Fraction


def canonical(forest):
    return tuple(sorted(forest, key=repr))


@lru_cache(None)
def vertex_count(tree):
    return 1 + sum(map(vertex_count, tree))


def parse_tree(text, ordinary=False):
    tree = ast.literal_eval(text)

    def valid(node):
        require((isinstance(node, tuple)), "A tree/forest must be a tuple.")
        if ordinary:
            require((node == canonical(node)), "Noncanonical ordinary tree.")
        for child in node:
            valid(child)

    valid(tree)
    return tree


@lru_cache(None)
def enumerate_trees(vertices, ordered):
    if vertices == 1:
        return ((),)
    results = set()

    def grow(remaining, children):
        if not remaining:
            results.add(children if ordered else canonical(children))
            return
        for count in range(1, remaining + 1):
            for child in enumerate_trees(count, ordered):
                grow(remaining - count, children + (child,))

    grow(vertices - 1, ())
    return tuple(sorted(results, key=repr))


def read_polynomial(records, weight):
    polynomial = {}
    for term in records:
        forest = parse_tree(term["forest"], ordinary=True)
        coefficient = Q(term["coefficient"])
        require((coefficient and forest not in polynomial), 'verify_certificate.py:70: check failed')
        require((sum(map(vertex_count, forest)) == weight), 'verify_certificate.py:71: check failed')
        polynomial[forest] = coefficient
    return polynomial


def combine(terms):
    answer = defaultdict(Q)
    for polynomial, multiplier in terms:
        for forest, coefficient in polynomial.items():
            answer[forest] += multiplier * coefficient
    return {forest: value for forest, value in answer.items() if value}


def multiply(left, right):
    answer = defaultdict(Q)
    for forest, coefficient in left.items():
        for other, value in right.items():
            answer[canonical(forest + other)] += coefficient * value
    return {forest: value for forest, value in answer.items() if value}


@lru_cache(None)
def all_words(multiset):
    return tuple(sorted(set(permutations(multiset)), key=repr))


@lru_cache(None)
def shuffle_multiplicities(left, right):
    """Enumerate interleavings by their chosen left-word positions."""
    counts = Counter()
    length = len(left) + len(right)
    for left_positions in combinations(range(length), len(left)):
        positions = set(left_positions)
        first, second = iter(left), iter(right)
        word = tuple(next(first) if i in positions else next(second)
                     for i in range(length))
        counts[word] += 1
    return dict(counts)


def contraction_polynomial(first, tail):
    """Expand b R_tail[-1] ... R_tail[0] Phi_first from the stated R_v."""
    terms = {((), first): Q(1)}
    for argument in tail:
        following = defaultdict(Q)
        for (forest, stage_tree), coefficient in terms.items():
            # (b Phi_argument) Phi_stage_tree
            following[(canonical(forest + (argument,)), stage_tree)] += coefficient
            # -(A Phi_argument) Phi_stage_tree
            following[(forest, canonical(stage_tree + (argument,)))] -= coefficient
            # +Phi_argument (A Phi_stage_tree)
            following[(forest, canonical(argument + (stage_tree,)))] += coefficient
        terms = {key: value for key, value in following.items() if value}
    result = defaultdict(Q)
    for (forest, stage_tree), coefficient in terms.items():
        result[canonical(forest + (stage_tree,))] += coefficient
    return {forest: value for forest, value in result.items() if value}


def matrix_rank(matrix):
    """Exact Gauss-Jordan rank, independent of the generator's solver."""
    if not matrix:
        return 0
    matrix = [list(map(Q, row)) for row in matrix]
    row = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(row, len(matrix))
                      if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        divisor = matrix[row][column]
        matrix[row] = [value / divisor for value in matrix[row]]
        for i in range(len(matrix)):
            if i != row and matrix[i][column]:
                multiplier = matrix[i][column]
                matrix[i] = [x - multiplier * y
                             for x, y in zip(matrix[i], matrix[row])]
        row += 1
        if row == len(matrix[0]):
            break
    return row


def equation_rows(words):
    rows, labels = [], []
    for first in sorted(set(words[0]), key=repr):
        rows.append([int(word[0] == first) for word in words])
        labels.append("P(" + str(first) + ")")
    for word in words:
        for cut in range(1, len(word)):
            left, right = word[:cut], word[cut:]
            counts = shuffle_multiplicities(left, right)
            rows.append([counts.get(other, 0) for other in words])
            labels.append("S(" + "".join(map(str, left)) + ","
                          + "".join(map(str, right)) + ")")
    return rows, labels


def leibniz_determinant(matrix):
    n = len(matrix)
    require((all(len(row) == n for row in matrix)), 'verify_certificate.py:172: check failed')
    result = 0
    for ordering in permutations(range(n)):
        sign = (-1) ** sum(ordering[i] > ordering[j]
                          for i in range(n) for j in range(i + 1, n))
        term = sign
        for row, column in enumerate(ordering):
            term *= matrix[row][column]
        result += term
    return result


@lru_cache(None)
def tree_factorial(tree):
    result = vertex_count(tree)
    for child in tree:
        result *= tree_factorial(child)
    return result


def at_exact_moments(polynomial):
    total = Q(0)
    for forest, coefficient in polynomial.items():
        for tree in forest:
            coefficient /= tree_factorial(tree)
        total += coefficient
    return total


@lru_cache(None)
def exact_flow_coefficient(tree):
    result, prefix_weight = Q(1), 0
    for child in tree:
        prefix_weight += vertex_count(child)
        result *= exact_flow_coefficient(child) / prefix_weight
    return result


def check_certificate(certificate_path, ranks_path):
    report = json.loads(certificate_path.read_text())
    require((report["status"] == "PASS" and report["order_requested"] == 6), 'verify_certificate.py:212: check failed')
    formulas = {}
    family_words = []
    for family in report["tensor_formulas"]:
        multiset = parse_tree(family["multiset"], ordinary=True)
        weight = sum(map(vertex_count, multiset))
        require((weight == family["weight"] <= 6), 'verify_certificate.py:218: check failed')
        expected = set(all_words(multiset))
        found = set()
        for entry in family["formulas"]:
            word = parse_tree(entry["word"])
            require((word not in formulas and word not in found), 'verify_certificate.py:223: check failed')
            require((canonical(word) == multiset), 'verify_certificate.py:224: check failed')
            polynomial = read_polynomial(entry["polynomial"], weight)
            formulas[word] = polynomial
            found.add(word)
        require((found == expected), 'verify_certificate.py:228: check failed')
        family_words.append(tuple(sorted(found, key=repr)))

    polarisation_checks = shuffle_checks = rank_checks = 0
    for words in family_words:
        length = len(words[0])
        if length == 1:
            require((formulas[words[0]] == {(words[0][0],): Q(1)}), 'verify_certificate.py:235: check failed')
            continue
        multiset = canonical(words[0])
        for first in set(multiset):
            tail = list(multiset)
            tail.remove(first)
            tails = all_words(tuple(tail))
            lhs = combine((formulas[(first,) + ordered_tail], Q(1))
                          for ordered_tail in tails)
            rhs = combine((contraction_polynomial(first, ordered_tail),
                           Q(1, factorial(length))) for ordered_tail in tails)
            require((lhs == rhs), ("Polarisation mismatch", multiset, first))
            polarisation_checks += 1
        for word in words:
            for cut in range(1, length):
                left, right = word[:cut], word[cut:]
                lhs = combine((formulas[shuffled], multiplicity)
                              for shuffled, multiplicity in
                              shuffle_multiplicities(left, right).items())
                rhs = multiply(formulas[left], formulas[right])
                require((lhs == rhs), ("Shuffle mismatch", left, right))
                shuffle_checks += 1
        rows, _ = equation_rows(words)
        require((matrix_rank(rows) == len(words)), 'verify_certificate.py:258: check failed')
        rank_checks += 1

    coefficients = {}
    records = report["coefficients"]
    for entry in records:
        tree = parse_tree(entry["tree"])
        degree = vertex_count(tree) - 1
        require((degree == entry["order"] and tree not in coefficients), 'verify_certificate.py:266: check failed')
        coefficients[tree] = read_polynomial(entry["polynomial"], degree)
    expected_trees = {tree for n in range(1, 8)
                      for tree in enumerate_trees(n, True)}
    require((set(coefficients) == expected_trees and len(coefficients) == 197), 'verify_certificate.py:270: check failed')

    for entry in records:
        tree = parse_tree(entry["tree"])
        if not tree:
            reconstructed = {(): Q(1)}
        else:
            choices = [list(coefficients[child].items()) for child in tree]
            summands = []
            for selection in product(*choices):
                arguments = tuple(forest for forest, _ in selection)
                multiplier = Q(1)
                for _, coefficient in selection:
                    multiplier *= coefficient
                summands.append((formulas[arguments], multiplier))
            reconstructed = combine(summands)
        require((reconstructed == coefficients[tree]), ("Root recursion", tree))
        actual = at_exact_moments(coefficients[tree])
        target = exact_flow_coefficient(tree)
        require((actual == target), 'verify_certificate.py:289: check failed')
        require((actual == Q(entry["exact_moment_value"])), 'verify_certificate.py:290: check failed')
        require((target == Q(entry["exact_flow_value"])), 'verify_certificate.py:291: check failed')

    ranks = json.loads(ranks_path.read_text())
    rank_summary = []
    for certificate in ranks:
        pattern = tuple(certificate["multiplicity_pattern"])
        words = tuple(sorted(set(permutations(pattern))))
        require((["".join(word) for word in words] == certificate["words"]), 'verify_certificate.py:298: check failed')
        rows, labels = equation_rows(words)
        available = dict(zip(labels, rows))
        asserted = certificate["matrix"]
        for label, row in zip(certificate["equation_labels"], asserted):
            require((available[label] == row), 'verify_certificate.py:303: check failed')
        computed_rank = matrix_rank(rows)
        require((computed_rank == certificate["rank"]), 'verify_certificate.py:305: check failed')
        if computed_rank == len(words):
            determinant = leibniz_determinant(asserted)
            require((determinant == Q(certificate["determinant"]) != 0), 'verify_certificate.py:308: check failed')
        else:
            require((pattern == tuple("aaabb") and computed_rank == 9), 'verify_certificate.py:310: check failed')
            require((len(words) == 10 and certificate["determinant"] is None), 'verify_certificate.py:311: check failed')
        rank_summary.append({"pattern": "".join(pattern),
                             "rank": computed_rank, "unknowns": len(words)})

    # Independently enumerate all ordinary-argument multisets of weight <= 6.
    ordinary = [tree for n in range(1, 7)
                for tree in enumerate_trees(n, False)]
    multiplicities = set()

    def classify(start, remaining, arguments):
        if arguments:
            shape = tuple(sorted(Counter(arguments).values(), reverse=True))
            if len(shape) > 1:
                multiplicities.add(shape)
        for index in range(start, len(ordinary)):
            weight = vertex_count(ordinary[index])
            if weight <= remaining:
                classify(index, remaining - weight,
                         arguments + (ordinary[index],))

    classify(0, 6, ())
    expected_shapes = {(1, 1), (2, 1), (1, 1, 1),
                       (3, 1), (2, 2), (4, 1)}
    require((multiplicities == expected_shapes), 'verify_certificate.py:334: check failed')
    return {"status": "PASS", "checker_imports": "Python standard library only",
            "ordered_tree_polynomials": len(coefficients),
            "tensor_formulas": len(formulas),
            "nontrivial_tensor_family_ranks": rank_checks,
            "polarisation_equations": polarisation_checks,
            "shuffle_equations": shuffle_checks,
            "root_recursions": len(coefficients),
            "exact_flow_values": len(coefficients),
            "small_rank_certificates": rank_summary,
            "multiplicity_classification_through_weight_six": "PASS",
            "order_seven_statement": "One undetermined direction; no counterexample asserted."}


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path,
                        default=directory / "order6_symbolic_certificate.json")
    parser.add_argument("--ranks", type=Path,
                        default=directory / "tensor_ranks.json")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = check_certificate(arguments.certificate, arguments.ranks)
    serialized = json.dumps(result, indent=2) + "\n"
    if arguments.output:
        arguments.output.write_text(serialized)
    print(serialized, end="")
