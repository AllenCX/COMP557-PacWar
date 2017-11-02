import _PyPacwar
import numpy as np
import random
import sys, os, json, time, copy
import operator
from multiprocessing import Pool, Process, Queue, Array

def generate_GENE(population, gene_length=50):
    ge = np.random.choice(4, gene_length).tolist()
    GE = []
    for i in xrange(population):
        GE.append(np.random.choice(4, gene_length).tolist())
    return GE

def run_pacWar(GE1, GE2, idx):
    (rounds, c1, c2) = _PyPacwar.battle(GE1[idx[0]], GE2[idx[1]]) 
    scores = GAPacWar.get_scores(rounds, c1, c2)
    ge1_scores[idx[0]] += scores[0]
    ge2_scores[idx[1]] += scores[1]
    #print(scores)
    return scores
class GAPacWar(object):
    def __init__(self, 
                 GE1,
                 GE2, 
                 mutation_prob=1e-2, 
                 num_iters=1000, 
                 elimination_step=200,
                 log_step=10,
                 save_step=100,
                 seed=7, 
                 run_id=None, 
                 save_dir=""): 
        self.run_id = run_id
        self.save_dir = save_dir
        self.GE1 = GE1 
        self.GE2 = GE2
        self.mutation_prob = mutation_prob
        self.num_iters = num_iters
        self.elimination_step = elimination_step
        self.log_step = log_step
        self.save_step = save_step
        self.seed = seed
        if not run_id:
            self.run_id = time.strftime("%y%m%d%H%M%S")
        else:
            self.run_id = run_id
        self.save_dir = save_dir
        np.random.seed(seed)
    @staticmethod
    def get_scores(rounds, c1, c2):
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


    def evaluate(self, GE1, GE2, parallel=True, verbose=0):
        '''
        GE1: genes for specie 1
        GE2: genes for specie 2
        '''
        #assert len(GE1) == len(GE2)
        #pop_size = len(GE1)
        #start_time = time.time()
        if parallel:
            assert verbose == 0
            global ge1_scores
            global ge2_scores
            ge1_scores = Array('d', [0] * len(GE1))
            ge2_scores = Array('d', [0] * len(GE2))
            
            p = Pool(4)
            for i in xrange(len(GE1)):
                for j in xrange(len(GE2)):
                    # bottleneck
                    try:
                        res = p.apply_async(run_pacWar, args=(GE1, GE2, (i, j),))
                        
                    #scores = res.get(timeout=1)
                    if verbose:
                        print("The duel lasted %d rounds, specie 1 remains: %d, specie 2 remains: %d" % (rounds, c1, c2))
                        print("The final scores for GE1[%d] vs GE2[%d] are: %d : %d" % (i, j, scores[0], scores[1]))
            try:
                p.close()
                p.join()
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                p.terminate()
                p.join()
            ge1_scores = list(ge1_scores)
            ge2_scores = list(ge2_scores)
            #print(time.time() - start_time)
        if not parallel:
            ge1_scores = [0.0] * len(GE1)
            ge2_scores = [0.0] * len(GE2)
            
            for i in xrange(len(GE1)):
                for j in xrange(len(GE2)):
                    # bottleneck
                    (rounds, c1, c2) = _PyPacwar.battle(GE1[i], GE2[j])
                    
                    scores = self.get_scores(rounds, c1, c2)
                    ge1_scores[i] += scores[0]
                    ge2_scores[j] += scores[1]
                    if verbose:
                        print("The duel lasted %d rounds, specie 1 remains: %d, specie 2 remains: %d" % (rounds, c1, c2))
                        print("The final scores for GE1[%d] vs GE2[%d] are: %d : %d" % (i, j, scores[0], scores[1]))
        ge1_scores = [i / len(GE2) for i in ge1_scores]
        ge2_scores = [i / len(GE1) for i in ge2_scores]
        return ge1_scores, ge2_scores

    def selection(self, GE1, GE2, s1, s2):
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

    def mate_and_mutate(self, GE1, GE2, s1, s2, mutation_prob):
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

    def train(self):
        GE1, GE2 = self.GE1, self.GE2
        mutation_prob = self.mutation_prob
        random_gene = [[3 for _ in xrange(50)]]

        #self.evaluate(GE1, random_gene, verbose=1)
        for i in xrange(1, self.num_iters+1):
            s1, s2 = self.evaluate(GE1, GE2)
            
            sys.stdout.write("At step %d, The average score for specie 1: %.2f, specie 2: %.2f\r" % 
                (i, sum(s1) / float(len(s1)), sum(s2) / float(len(s2))))
            '''
            if i % self.log_step == 0:
                print("At step %d, The average score for specie 1: %.2f, specie 2: %.2f" % 
                    (i, sum(s1) / float(len(s1)), sum(s2) / float(len(s2))))
            '''
            if i % self.save_step == 0:
                self.save_checkpoint(GE1, GE2, s1, s2, self.run_id, i, self.save_dir)
                print("ckpt saved!")
            '''
            if i % self.elimination_step == 0 and i != 0:
                total_s1, total_s2 = sum(s1), sum(s2)
                if total_s1 > total_s2:
                    GE2 = copy.deepcopy(GE1)
                    print("Specie 2 eliminated!")
                else:
                    GE1 = copy.deepcopy(GE2)
                    print("Specie 1 eliminated!")
            '''
            s1, s2, GE1, GE2 = self.selection(GE1, GE2, s1, s2)
            GE1, GE2 = self.mate_and_mutate(GE1, GE2, s1, s2, mutation_prob)
            
        #(rounds, c1, c2) = _PyPacwar.battle(GE1[0], random_gene[0])
        #print("The duel lasted %d rounds, specie 1 remains: %d, random specie remains: %d" % (rounds, c1, c2))
        self.evaluate(GE1, random_gene, verbose=1)
        self.evaluate(GE2, random_gene, verbose=1)

        best_gene, best_score = self.find_the_best(GE1, GE2)
        print(best_gene, best_score)
        return GE1, GE2

    def find_the_best(self, GE1, GE2, verbose=0):
        ge1_scores = [0] * len(GE1)
        ge2_scores = [0] * len(GE2)
        for i in xrange(len(GE1)):
            for j in xrange(len(GE2)):
                (rounds, c1, c2) = _PyPacwar.battle(GE1[i], GE2[j])
                scores = self.get_scores(rounds, c1, c2)
                #ge1_scores[i] += scores[0]
                #ge2_scores[j] += scores[1]
                if verbose:
                    print("The duel lasted %d rounds, specie 1 remains: %d, specie 2 remains: %d" % (rounds, c1, c2))
                    print("The final scores for GE1[%d] vs GE2[%d] are: %d : %d" % (i, j, scores[0], scores[1]))
        best_index1, best_score_1 = max(enumerate(ge1_scores), key=operator.itemgetter(1))
        best_index2, best_score_2 = max(enumerate(ge2_scores), key=operator.itemgetter(1))
        #return GE1[best_index1], GE2[best_index2], best_score_1, best_score_2
        if best_score_1 > best_score_2:
            return GE1[best_index1], best_score_1
        else:
            return GE2[best_index2], best_score_2

    def test(self, GE1, GE2):
        pass

    def save_checkpoint(self, GE1, GE2, s1, s2, run_id, step, save_dir):
        '''
        save checkpoint and find the best single gene in this two species
        '''
        GE1_dict, GE2_dict = {}, {}
        for i in xrange(len(GE1)):
            gene_str = "".join([str(k) for k in GE1[i]])
            GE1_dict[gene_str] = s1[i]
        for i in xrange(len(GE2)):
            gene_str = "".join([str(k) for k in GE2[i]])
            GE2_dict[gene_str] = s2[i]
        GE1_json_str, GE2_json_str = json.dumps(GE1_dict), json.dumps(GE2_dict)
        with open("%s/%s_GE1.json" % (save_dir, str(run_id)), "w") as f:
            json.dump(GE1_dict, f)
        with open("%s/%s_GE2.json" % (save_dir, str(run_id)), "w") as f:
            json.dump(GE2_dict, f)
        best_gene, best_score = self.find_the_best(GE1, GE2)

if __name__ == '__main__':
    GE1 = [[3 for _ in xrange(50)] for _ in xrange(50)]
    GE2 = [[3 for _ in xrange(50)] for _ in xrange(50)]
    gaParwar = GAPacWar(GE1, GE2, mutation_prob=0.02, num_iters=100, save_step=10, run_id="50pop", save_dir="ckpt/")
    gaParwar.train()