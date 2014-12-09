#include "globheads.h"
#include "protos.h"
#include <assert.h>
#include <stdio.h>

typedef struct {
    int len; // amount of nodes
    int *nodes;
    int size; //size of nodes
    int last; //last node
    int ref; //ref count
    int processed;
} Clique;

typedef struct {
    int len; // amount of cliques
    Clique **cliques;
    int size; // size of cliques
    Clique **current;
} Cliques;

void nodes_equal(Clique *c, int *nodes, int len)
{
    static int try = 0;
    try++;
    int i;
    printf("try %d, len(c) %d, len %d\n", try, c->len, len);
    assert(len == c->len);
    for (i = 0; i < len; ++i)
    {
        if (nodes[i] != c->nodes[i])
        {
            printf("%d != %d at %d try %d\n", nodes[i], c->nodes[i], i, try);
            clique_print(c);
            assert(nodes[i] == c->nodes[i]);
        }
    }
}

int main()
{
    int nodes[5] = {2, 6, 33, 34, 63};
    Clique *c = clique_new_from_nodes(nodes, 5);
    nodes_equal(c, nodes, 5);

    Clique *c2 = clique_new_from_clique(c);
    nodes_equal(c2, nodes, 5);

    int nodes2[6] = {2, 6, 15, 33, 34, 63};
    clique_insert(c2, 15);
    nodes_equal(c2, nodes2, 6);

    Clique *c3 = clique_intersection(c, c2);
    nodes_equal(c3, nodes, 5);
    clique_free(c3);

    c3 = clique_union(c, c2);
    nodes_equal(c3, nodes2, 6);
    clique_free(c3);

    c3 = clique_exclude(c, c2);
    nodes_equal(c3, NULL, 0);
    clique_free(c3);

    int nodes3[1] = {15};
    c3 = clique_exclude(c2, c);
    nodes_equal(c3, nodes3, 1);
    clique_free(c3);

    c3 = clique_new_from_clique(c2);
    clique_remove(c3, 15);
    nodes_equal(c3, nodes, 5);

    int nodes4[6] = {1, 2, 6, 33, 34, 63};
    clique_insert(c3, 1);
    nodes_equal(c3, nodes4, 6);

    int nodes5[7] = {1, 2, 6, 33, 34, 45, 63};
    clique_insert(c3, 45);
    nodes_equal(c3, nodes5, 7);

    int nodes6[8] = {1, 2, 6, 33, 34, 45, 63, 75};
    clique_insert(c3, 75);
    nodes_equal(c3, nodes6, 8);

    int nodes7[9] = {1, 2, 6, 33, 34, 45, 63, 75, 80};
    clique_append(c3, 80);
    nodes_equal(c3, nodes7, 9);

    int nodes8[10] = {1, 2, 6, 15, 33, 34, 45, 63, 75, 80};
    clique_remove(c3, 33);
    Clique *c4 = clique_union(c3, c2);
    nodes_equal(c4, nodes8, 10);

    clique_remove(c4, 12);
    nodes_equal(c4, nodes8, 10);
    clique_remove(c4, 0);
    nodes_equal(c4, nodes8, 10);
    clique_remove(c4, 112);
    nodes_equal(c4, nodes8, 10);
    clique_free(c4);

    clique_free(c3);

    nodes_equal(c, nodes, 5);
    nodes_equal(c2, nodes2, 6);

    clique_clear(c2);
    nodes_equal(c2, NULL, 0);

    clique_free(c2);
    clique_free(c);
    return 0;
}
