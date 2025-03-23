import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    probability = 1.0
    for person in people:
        # Define the number of genes
        gene = 0
        if person in one_gene:
            gene = 1
        elif person in two_genes:    
            gene = 2
        else:
            gene = 0   

        # Get the names of the parents
        mother = people[person]['mother']
        father = people[person]['father']

        probability_gene = 0.0
        # Calculates the probability of the gene
        if not mother and not father: # If the person has no parents
            probability_gene = PROBS['gene'][gene] # unconditional probability

        else:
            # Create parents info objects  
            parents_info = {
                # For each parent
                parent: {
                    # Define the number of genes
                    'gene': 1 if parent in one_gene else 2 if parent in two_genes else 0, 

                    # Define the probability of the gene to be passed (0.0 for now)
                    'gene_inheritance': 0.0, 
                } for parent in [mother, father]
            }
            
            # Calculate the probability of pass the gene from parents
            for parent in parents_info:

                # Get the number of genes
                parent_gene = parents_info[parent]['gene']

                # Get the probability of the gene to be passed
                if parent_gene == 0:
                    parents_info[parent]['gene_inheritance'] = PROBS['mutation'] # Normal prob
                elif parent_gene == 1:
                    # The probability of the gene to be passed is 50% of one parent and other 50% of the other parent
                    # First 50% with normal probability + 50% with mutation probability 
                    gene_inheritance_value = 0.5 * PROBS['mutation'] + 0.5 * (1 - PROBS['mutation'])
                    parents_info[parent]['gene_inheritance'] = gene_inheritance_value
                else:
                    # The probability of two genes to be passed is 1 / mutation probability
                    # 1 are both = 0.5 one parent + 0.5 other parent and all that divided by mutation probability 
                    parents_info[parent]['gene_inheritance'] = 1 - PROBS['mutation']

            # Calculate probability of gene keeping in mind the inheritance probs
            if gene == 0:
                probability_gene = (1 - parents_info[mother]['gene_inheritance']) * \
                (1 - parents_info[father]['gene_inheritance'])
            elif gene == 1:
                probability_gene = ((1 - parents_info[mother]['gene_inheritance']) * parents_info[father]['gene_inheritance'] +
                                    (1 - parents_info[father]['gene_inheritance']) * parents_info[mother]['gene_inheritance'])
            else:
                probability_gene = parents_info[mother]['gene_inheritance'] * parents_info[father]['gene_inheritance']
            
            # Calculate trait probability
            probability_trait = PROBS['trait'][gene][True if person in have_trait else False]

            probability *= probability_gene * probability_trait
            print(probability)
            return probability



def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    raise NotImplementedError


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
