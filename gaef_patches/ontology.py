"""GO ancestry for evaluation: is_a and part_of, never regulates.

Replace the former import-time GAEF monkeypatch with explicit IC methods.
Reuse the bundled DeepGO implementation with explicit relation semantics and
defined behavior for a corpus whose information contents are all zero.
"""
from collections import Counter
import math
from deepgo.utils import Ontology as BaseOntology, FUNC_DICT, NAMESPACES


class Ontology(BaseOntology):
    def __init__(self, filename='data/go-basic.obo', with_rels=True):
        self.filename = str(filename)
        super().__init__(filename, with_rels=False)
        self.aliases = {key: obj['id'] for key, obj in self.ont.items()
                        if not obj['is_obsolete']}
        self.ont = {key: obj for key, obj in self.ont.items()
                    if key == obj['id'] and not obj['is_obsolete']}
        for obj in self.ont.values():
            obj['is_a'] = sorted({self.aliases[parent] for parent in obj['is_a']
                                 if self.aliases.get(parent) in self.ont})
            obj['children'] = set()
        for term, obj in self.ont.items():
            for parent in obj['is_a']:
                self.ont[parent]['children'].add(term)
        if with_rels:
            term = None
            with open(filename) as fh:
                for raw in fh:
                    line = raw.strip()
                    if line.startswith('['):
                        term = None
                    elif line.startswith('id: GO:'):
                        term = line.split()[1]
                    elif term in self.ont and line.startswith('relationship: part_of GO:'):
                        parent = self.canonical_id(line.split()[2])
                        if parent in self.ont:
                            self.ont[term]['is_a'].append(parent)
                            self.ont[parent]['children'].add(term)
        self.relations = ('is_a', 'part_of') if with_rels else ('is_a',)

    def canonical_id(self, term):
        primary = self.aliases.get(term)
        return primary if primary in self.ont else None

    def normalize_terms(self, terms):
        return {primary for term in terms if (primary := self.canonical_id(term)) is not None}

    def normalize_scores(self, predictions, namespace=None, audit=None):
        result = {}
        counts = Counter()
        for protein, scores in predictions.items():
            normalized = {}
            for term, score in scores.items():
                if not math.isfinite(score) or not 0 <= score <= 1:
                    raise ValueError(f'Invalid prediction score: {protein}, {term}, {score}')
                primary = self.canonical_id(term)
                if primary is None:
                    counts['unknown_or_obsolete_pairs'] += 1
                    continue
                if namespace is not None and self.ont[primary]['namespace'] != namespace:
                    counts['other_namespace_pairs'] += 1
                    continue
                counts['alternative_id_pairs'] += int(primary != term)
                counts['merged_duplicate_pairs'] += int(primary in normalized)
                normalized[primary] = max(score, normalized.get(primary, 0.0))
            result[protein] = normalized
        if audit is not None:
            audit.update(counts)
        return result

    def normalize_constraints(self, mapping):
        result = {}
        for term, values in mapping.items():
            primary = self.canonical_id(term)
            if primary is not None:
                result.setdefault(primary, set()).update(values)
        return result

    def propagate_scores(self, predictions, audit=None):
        """Max-propagate canonical scores before comparing repair methods."""
        result = {}
        added = raised = 0
        for protein, scores in predictions.items():
            closed = dict(scores)
            for term, score in scores.items():
                if score <= 0:
                    continue
                for ancestor in self.get_ancestors(term):
                    closed[ancestor] = max(score, closed.get(ancestor, 0.0))
            added += sum(g not in scores for g in closed)
            raised += sum(g in scores and value > scores[g] for g, value in closed.items())
            result[protein] = closed
        if audit is not None:
            audit.update(hierarchy_added_pairs=added, hierarchy_raised_pairs=raised)
        return result

    def get_ancestors(self, term):
        return super().get_ancestors(self.canonical_id(term))

    def get_parents(self, term):
        return super().get_parents(self.canonical_id(term))

    def get_term_set(self, term):
        return super().get_term_set(self.canonical_id(term))

    def get_ic(self, term):
        return super().get_ic(self.canonical_id(term))

    def has_term(self, term):
        return self.canonical_id(term) is not None

    def get_term(self, term):
        return self.ont.get(self.canonical_id(term))

    def get_namespace(self, term):
        return self.ont[self.canonical_id(term)]['namespace']

    def calculate_ic(self, annots):
        counts = Counter()
        for terms in annots:
            closed = set()
            for term in terms:
                closed.update(self.get_ancestors(term))
            counts.update(closed)
        self.calculate_ic_from_counts(counts)

    def calculate_ic_from_counts(self, counts):
        """Use a fixed ancestor-closed reference corpus for every comparator."""
        self.ic = {}
        self.ic_norm = 0.0
        for term, n in counts.items():
            parents = self.get_parents(term)
            parent_n = min((counts.get(p, 0) for p in parents), default=n)
            if parent_n < n:
                raise ValueError(f'IC corpus is not ancestor-closed at {term}')
            value = math.log2(parent_n / n)
            self.ic[term] = value
            self.ic_norm = max(self.ic_norm, value)

    def get_norm_ic(self, go_id):
        return self.get_ic(go_id) / self.ic_norm if self.ic_norm else 0.0
