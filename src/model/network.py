from model.disease import Disease
from model.drug import Drug
from model.gene import Gene
from model.variant import Variant


class Network:
    def __init__(self):
        self.drugs = {}
        self.genes = {}
        self.variants = {}
        self.diseases = {}
        self.drug_indicates_disease = []
        self.drug_contraindicates_disease = []
        self.drug_targets_gene = []
        self.gene_associates_with_disease = []
        self.variant_associates_with_disease = []

    def add_drug(self, drug: Drug):
        matches = [self.drugs[x] for x in drug.ids if x in self.drugs]
        for match in matches:
            drug.merge(match)
        for x in drug.ids:
            self.drugs[x] = drug

    def add_gene(self, gene: Gene):
        matches = [self.genes[x] for x in gene.ids if x in self.genes]
        for match in matches:
            gene.merge(match)
        for x in gene.ids:
            self.genes[x] = gene

    def add_disease(self, disease: Disease):
        matches = [self.diseases[x] for x in disease.ids if x in self.diseases]
        for match in matches:
            disease.merge(match)
        for x in disease.ids:
            self.diseases[x] = disease

    def add_variant(self, variant: Variant):
        matches = [self.variants[x] for x in variant.ids if x in self.variants]
        for match in matches:
            variant.merge(match)
        for x in variant.ids:
            self.variants[x] = variant

    def get_drug_by_id(self, _id: str) -> Drug:
        return self.drugs[_id]

    def get_disease_by_id(self, _id: str) -> Disease:
        return self.diseases[_id]

    def get_gene_by_id(self, _id: str) -> Gene:
        return self.genes[_id]

    def get_variant_by_id(self, _id: str) -> Variant:
        return self.variants[_id]

    def prune(self):
        """
        Remove nodes that are not connected to drugs and therefore of no interest.
        """
        removed_gene_ids = set()
        targeted_genes_id = {x[1] for x in self.drug_targets_gene}
        for gene in set(self.genes.values()):
            if targeted_genes_id.isdisjoint(gene.ids):
                for gene_id in gene.ids:
                    del self.genes[gene_id]
                    removed_gene_ids.add(gene_id)
