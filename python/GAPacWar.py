import _PyPacwar
import numpy as np
import random
import sys, os, json
import operator

def generate_GENE(population, gene_length=50):
    ge = np.random.choice(4, gene_length).tolist()
    GE = []
    for i in xrange(population):
        GE.append(np.random.choice(4, gene_length).tolist())
    return GE


def get_scores(rounds, c1, c2):
    '''
    rounds: the real rounds which the duel last
    c1: our species remaining number
    c2: oppenents' species remaining number
    '''
    if c1 == c2:
        return [10, 10]

    species = [c1, c2]
    scores = [0, 0]
    winner = 0 if c1 > c2 else 1
    loser = abs(1 - winner)

    if species[loser] == 0:
        if rounds < 100:
            scores[winner], scores[loser] = 20, 0 
        elif 100 <= rounds <= 199:
            scores[winner], scores[loser] = 19, 1
        elif 200 <= rounds <= 299:
            scores[winner], scores[loser] = 18, 2
        elif 300 <= rounds <= 500:
            scores[winner], scores[loser] = 17, 3
    else:
        if species[winner] >= 10 * species[loser]:
            scores[winner], scores[loser] = 13, 7
        elif 3 * species[loser] <= species[winner] < 10 * species[loser]:
            scores[winner], scores[loser] = 12, 8
        elif 1.5 * scores[loser] <= species[winner] < 3 * species[loser]:
            scores[winner], scores[loser] = 11, 9
        else:
            scores[winner], scores[loser] = 10, 10
    return scores


def evaluation(GE1, GE2, verbose=0):
    '''
    GE1: genes for specie 1
    GE2: genes for specie 2
    '''
    #assert len(GE1) == len(GE2)
    #pop_size = len(GE1)
    ge1_scores = [0] * len(GE1)
    ge2_scores = [0] * len(GE2)
    for i in xrange(len(GE1)):
        for j in xrange(len(GE2)):
            (rounds, c1, c2) = _PyPacwar.battle(GE1[i], GE2[j])
            scores = get_scores(rounds, c1, c2)
            ge1_scores[i] += scores[0]
            ge2_scores[j] += scores[1]
            if verbose:
                print("The duel lasted %d rounds, specie 1 remains: %d, specie 2 remains: %d" % (rounds, c1, c2))
                print("The final scores for GE1[%d] vs GE2[%d] are: %d : %d" % (i, j, scores[0], scores[1]))
            
    # normalization
    return ge1_scores, ge2_scores

def selection(s1, s2, GE1, GE2):
    '''
    s1: scores for species 1
    s2: scores for species 2
    '''
    assert len(GE1) == len(GE2)
    pop_size = len(GE1)
    
    sum_ge1_scores = float(sum(s1)) + 1e-8
    sum_ge2_scores = float(sum(s2)) + 1e-8
    s1 = [i / sum_ge1_scores for i in s1]
    s2 = [i / sum_ge2_scores for i in s2]

    new_gene1_idx = np.random.choice(pop_size, size=pop_size, p=s1)
    new_gene1 = [GE1[i] for i in new_gene1_idx]
    new_s1 = [s1[i] for i in new_gene1_idx]

    new_gene2_idx = np.random.choice(pop_size, size=pop_size, p=s2)
    new_gene2 = [GE2[i] for i in new_gene2_idx]
    new_s2 = [s2[i] for i in new_gene2_idx]

    return new_s1, new_s2, new_gene1, new_gene2
    
def mate_and_mutate(s1, s2, GE1, GE2, mutation_prob):
    '''
    s1: a list scores of GE1
    s2: a list scores of GE2
    '''
    pop_size = len(GE1)
    scores_gene1 = zip(s1, GE1)
    scores_gene2 = zip(s2, GE2)

    # sort the scores by descend order
    scores_gene1 = sorted(scores_gene1, key=lambda x: x[0], reverse=True)
    scores_gene2 = sorted(scores_gene2, key=lambda x: x[0], reverse=True)
    
    # duplicate the gene with the highest score, remove the gene with lowest score
    scores_gene1.pop()
    scores_gene1.append(scores_gene1[0])
    scores_gene2.pop()
    scores_gene2.append(scores_gene2[0])
    
    prob1, gene1 = zip(*scores_gene1)
    prob2, gene2 = zip(*scores_gene2)

    # normalize
    sum_prob1 = float(sum(prob1))
    sum_prob2 = float(sum(prob2))
    prob1 = [i / float(sum(prob1)) for i in prob1]
    prob2 = [i / float(sum(prob2)) for i in prob2]


    def cross_over(gene1, gene2, pivot):
        '''
        gene1: a list represent gene
        gene2: a list represent gene
        pivot: cross over point, two genes will exchange the parts after the pivot
        '''
        gene1, gene2 = gene1[:pivot] + gene2[pivot:], gene2[:pivot] + gene1[pivot:]
        return gene1, gene2

    def mutate(gene, mutation_prob):
        for i in xrange(len(gene)):
            if np.random.rand() < mutation_prob:
                gene[i] = np.random.choice(4) # 4 is the number of different kinds of acids
        return gene
    # a list of shuffled genes, for each two (i, i + 1) of the elements, are the parents
    # which will mate and produce next generation
    parents1_idx, parents2_idx = range(pop_size), range(pop_size)
    np.random.shuffle(parents1_idx)
    np.random.shuffle(parents2_idx)

    #parents1_idx = np.random.choice(pop_size, pop_size, replace=False)
    #parents2_idx = np.random.choice(pop_size, pop_size, replace=False)

    for i in xrange(0, len(parents1_idx) - 1, 2):
        cross_over_pivot = np.random.choice(len(GE1[i]))
        GE1[i], GE1[i + 1] = cross_over(GE1[i], GE1[i+1], cross_over_pivot)
        GE1[i], GE1[i + 1] = mutate(GE1[i], mutation_prob), mutate(GE1[i+1], mutation_prob)

        cross_over_pivot = np.random.choice(len(GE2[i]))
        GE2[i], GE2[i + 1] = cross_over(GE2[i], GE2[i+1], cross_over_pivot)
        GE2[i], GE2[i + 1] = mutate(GE2[i], mutation_prob), mutate(GE2[i+1], mutation_prob)
    
    # Score is not necessary anymore here
    return GE1, GE2

def train(GE1, 
          GE2, 
          mutation_prob, 
          num_iters=2000, 
          log_step=10, 
          save_step=100, 
          run_id="0id",
          seed=7):
    '''
    GE1, GE2: initial genes for two species
    num_iters: number of iterations
    seed: random seed
    '''
    np.random.seed(7)
    print("Start training...")
    #random_gene = generate_GENE(1)
    random_gene = [[3 for _ in xrange(50)]]
    evaluation(GE1, random_gene, verbose=1)
    for i in xrange(num_iters):
        s1, s2 = evaluation(GE1, GE2)
        if i % log_step == 0:
           print("At step %d, The average score for specie 1: %.2f, specie 2: %.2f" % 
                (i, sum(s1) / float(len(s1)), sum(s2) / float(len(s2))))
        if i % save_step == 0:
            save_checkpoint(GE1, GE2, s1, s2, run_id, i)
        s1, s2, GE1, GE2 = selection(s1, s2, GE1, GE2)
        GE1, GE2 = mate_and_mutate(s1, s2, GE1, GE2, mutation_prob)
        

    
    #(rounds, c1, c2) = _PyPacwar.battle(GE1[0], random_gene[0])
    #print("The duel lasted %d rounds, specie 1 remains: %d, random specie remains: %d" % (rounds, c1, c2))
    evaluation(GE1, random_gene, verbose=1)
    evaluation(GE2, random_gene, verbose=1)
    return GE1, GE2

def find_the_best(GE1, GE2, verbose=0):
    ge1_scores = [0] * len(GE1)
    ge2_scores = [0] * len(GE2)
    for i in xrange(len(GE1)):
        for j in xrange(len(GE2)):
            (rounds, c1, c2) = _PyPacwar.battle(GE1[i], GE2[j])
            scores = get_scores(rounds, c1, c2)
            ge1_scores[i] += scores[0]
            ge2_scores[j] += scores[1]
            if verbose:
                print("The duel lasted %d rounds, specie 1 remains: %d, specie 2 remains: %d" % (rounds, c1, c2))
                print("The final scores for GE1[%d] vs GE2[%d] are: %d : %d" % (i, j, scores[0], scores[1]))
    best_index1, best_score_1 = max(enumerate(ge1_scores), key=operator.itemgetter(1))
    best_index2, best_score_2 = max(enumerate(ge2_scores), key=operator.itemgetter(1))
    return GE1[best_index1], GE2[best_index2], best_score_1, best_score_2

def save_checkpoint(GE1, GE2, s1, s2, run_id="0id", step=0):
    pass

def main():
    #np.random.seed(7)
    # generate two sets of genes
    #GE1 = generate_GENE(4)
    #GE2 = generate_GENE(4)
    GE1 = [[3 for _ in xrange(50)] for _ in xrange(10)]
    GE2 = [[3 for _ in xrange(50)] for _ in xrange(10)]
    GE1, GE2 = train(GE1, GE2, 1e-2)
    best_gene1, best_gene2, best_score_1, best_score_2 = find_the_best(GE1, GE2)
    print("".join([str(i) for i in best_gene1]), best_score_1)
    print("".join([str(i) for i in best_gene2]), best_score_2)
    #s1, s2 = evaluation(GE1, GE2)
    #s1, s2, GE1, GE2 = selection(s1, s2, GE1, GE2)
   
    #GE1, GE2 = mate_and_mutate(s1, s2, GE1, GE2, 1e-2)
    
    #(rounds,c1,c2) = _PyPacwar.battle(GE1, GE2)
    #print("The duel lasted %d rounds, specie 1 remains: %d, specie 2 remains: %d" % (rounds, c1, c2))
    #scores = get_scores(rounds, c1, c2)
    #print("The final scores are: %d : %d" % (scores[0], scores[1]))
    

if __name__ == "__main__":
    main()