#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "globheads.h"
#include "protos.h"
#include "profile.h"

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

int KeyComp( const void *vfst, const void *vsnd )
{
  KeyType *fst = (KeyType *)vfst, *snd = (KeyType *)vsnd;
  if( fst->key == snd->key ) {
    if( fst->var < snd->var )
      return -1;
    return 1;
  }
  if( fst->key < snd->key )
    return -1;
  return 1;
}

typedef struct {
    int len; // amount of nodes
    int *nodes;
    int size; //size of nodes
    int last; //last node
    int ref; //ref count
    int processed;
    void *neighbours;
    void *neighbour_blocks;
} Clique;

typedef struct {
    int len; // amount of cliques
    Clique **cliques;
    int size; // size of cliques
    Clique **current;
} Cliques;

Clique *clique_new()
{
    Clique *c = (Clique *) Malloc(sizeof(Clique), "clique_new");
    c->len = 0;
    c->size = 0;
    c->nodes = NULL;
    c->last = -1;
    c->ref = 1;
    c->processed = 0;
    c->neighbours = NULL;
    c->neighbour_blocks = NULL;
    return c;
}

Clique *clique_new_from_nodes(int *nodes, int len)
{
    Clique *c = (Clique *) Malloc(sizeof(Clique), "clique_new_from_nodes");
    c->len = len;
    c->size = len;
    c->nodes = (int *) Malloc(c->size*sizeof(int), "clique_new_from_nodes");
    memcpy(c->nodes, nodes, len*sizeof(int));
    c->last = c->nodes[len-1];
    c->ref = 1;
    c->processed = 0;
    c->neighbours = NULL;
    c->neighbour_blocks = NULL;
    return c;
}

void clique_unref(Clique *c);

void clique_free(Clique *c)
{
    if (c->nodes)
        free(c->nodes);
    if (c->neighbours)
        clique_unref(c->neighbours);
    if (c->neighbour_blocks)
        clique_unref(c->neighbour_blocks);
    free(c);
}

void clique_ref(Clique *c)
{
    c->ref++;
}

void clique_unref(Clique *c)
{
    if (!--c->ref)
        clique_free(c);
}

void clique_remove_neighbours(Clique *c)
{
    clique_unref((Clique *)c->neighbours);
    c->neighbours = NULL;
}

void clique_add_neighbours(Clique *c, Clique *n)
{
    if (c->neighbours)
        clique_remove_neighbours(c);

    c->neighbours = (void *)n;
    clique_ref(n);
}

void clique_remove_neighbour_blocks(Clique *c)
{
    clique_unref((Clique *)c->neighbour_blocks);
    c->neighbour_blocks = NULL;
}

void clique_add_neighbour_blocks(Clique *c, Clique *n)
{
    if (c->neighbour_blocks)
        clique_remove_neighbour_blocks(c);

    c->neighbour_blocks = (void *)n;
    clique_ref(n);
}

void clique_copy(Clique *dest, Clique *src)
{
    if (dest->nodes)
        free(dest->nodes);

    if (dest->neighbours)
        clique_remove_neighbours(dest);

    if (dest->neighbour_blocks)
        clique_remove_neighbour_blocks(dest);

    dest->len = src->len;
    dest->size = src->size;
    dest->last = src->last;
    dest->nodes = (int *) Malloc(dest->size*sizeof(int), "clique_copy");
    memcpy(dest->nodes, src->nodes, dest->size*sizeof(int));

    if (src->neighbours)
        clique_add_neighbours(dest, (Clique *)src->neighbours);

    if (src->neighbour_blocks)
        clique_add_neighbour_blocks(dest, (Clique *)src->neighbour_blocks);
}

Clique *clique_new_from_clique(Clique *src)
{
    Clique *dest = clique_new();
    clique_copy(dest, src);
    return dest;
}

void clique_clear(Clique *c)
{
    if (c->nodes)
        free(c->nodes);

    c->nodes = NULL;
    c->len = 0;
    c->size = 0;
    c->last = -1;
}
void clique_append(Clique *c, int v)
{
    if (c->last == v)
        return;

    if (c->size <= c->len)
    {
        c->size = max(1, c->size * 2);
        int *new_clique = (int *) Malloc(c->size*sizeof(int), "clique_append");
        if (c->nodes != NULL)
        {
            memcpy(new_clique, c->nodes, c->len*sizeof(int));
            free(c->nodes);
        }
        c->nodes = new_clique;
    }
    c->nodes[c->len++] = v;
    c->last = v;
}

void clique_append_clique(Clique *x, Clique *y)
{
    if (!y->len)
        return;

    if (x->size < x->len + y->len)
    {
        x->size = x->len + y->len;
        int *new_clique = (int *) Malloc(x->size*sizeof(int), "clique_append");
        if (x->nodes != NULL)
        {
            memcpy(new_clique, x->nodes, x->len*sizeof(int));
            free(x->nodes);
        }
        x->nodes = new_clique;
    }
    memcpy(x->nodes+x->len, y->nodes, y->len*sizeof(int));
    x->len += y->len;
    x->last = y->last;
}

int clique_bisect(Clique *x, int v)
{
    // returns index greater than or equal to v
    if (!x->len)
        return 0;

    int a = 0;
    int b = x->len - 1;

    if (x->nodes[a] >= v)
        return 0;

    if (x->nodes[b] < v)
        return x->len;

    int c = 0;
    while (b - a > 1)
    {
        c = a + (b - a) / 2;
        if (x->nodes[c] < v)
            a = c;
        else if (x->nodes[c] > v)
            b = c;
        else
            break;
    }

    if (v <= x->nodes[a])
        return a;

    if (v <= x->nodes[c])
        return c;

    return c+1;
}

void clique_print(Clique *c)
{
    if (!c->len)
    {
        printf("()\n");
        return;
    }

    printf("(");
    int iter;
    printf("%d", c->nodes[0]);
    for (iter = 1; iter < c->len; ++iter)
        printf(",%d", c->nodes[iter]);
    printf(")\n");
}

void clique_remove(Clique *x, int v)
{
    if (!x->len)
        return;

    int pos = clique_bisect(x, v);
    if (pos < x->len - 1 && x->nodes[pos] == v)
    {
        int size = (x->len - pos - 1) * sizeof(int);
        int *tmp = Malloc(size, "clique_remove");
        memcpy(tmp, x->nodes + pos + 1, size);
        memcpy(x->nodes + pos, tmp, size);
        x->len--;
        free(tmp);
    }
    else if (pos == x->len - 1 && x->nodes[pos] == v)
        if (--x->len)
            x->last = x->nodes[x->len-1];
}

void clique_insert(Clique *x, int v)
{
    if (!x->len)
    {
        clique_append(x, v);
        return;
    }

    int pos = clique_bisect(x, v);
    if (pos < x->len && x->nodes[pos] != v)
    {
        // Copy one less because we append that later
        int size = (x->len - pos - 1) * sizeof(int);
        int *tmp = Malloc(size, "clique_insert");
        memcpy(tmp, x->nodes + pos, size);
        memcpy(x->nodes + pos + 1, tmp, size);
        x->nodes[pos] = v;
        free(tmp);

        int tmp_last = x->last;
        x->last = -1;
        clique_append(x, tmp_last);
    }
    else if (pos == x->len)
        clique_append(x, v);
}

Clique *clique_union(Clique *x, Clique *y)
{
    Clique *out = clique_new();
    if (!x->len)
    {
        clique_copy(out, y);
        return out;
    }
    if (!y->len)
    {
        clique_copy(out, x);
        return out;
    }
    int i = 0;
    int j = 0;
    while (i < x->len && j < y->len)
    {
        if (x->nodes[i] < y->nodes[j])
            clique_append(out, x->nodes[i++]);
        else 
            clique_append(out, y->nodes[j++]);
    }
    while (i < x->len)
        clique_append(out, x->nodes[i++]);
    while (j < y->len)
        clique_append(out, y->nodes[j++]);
    return out;
}

Clique *clique_intersection(Clique *x, Clique *y)
{
    Clique *out = clique_new();
    if (!x->len || !y->len)
        return out;

    int i = 0;
    int j = 0;
    while (i < x->len && j < y->len)
    {
        if (y->nodes[j] == x->nodes[i])
        {
            clique_append(out, x->nodes[i]);
            i++;
            j++;
        }
        else if (x->nodes[i] < y->nodes[j])
            i++;
        else
            j++;
    }
    return out;
}

Clique *clique_exclude(Clique *x, Clique *y)
{
    Clique *out = clique_new();
    if (!x->len || !y->len)
    {
        clique_copy(out, x);
        return out;
    }
    int i = 0;
    int j = 0;
    while (i < x->len && j < y->len)
    {
        if (x->nodes[i] < y->nodes[j])
            clique_append(out, x->nodes[i++]);
        else if (x->nodes[i] > y->nodes[j])
            j++;
        else
            i++;
    }
    while (i < x->len)
    {
        if (x->nodes[i] != y->last)
            clique_append(out, x->nodes[i]);
        i++;
    }
    return out;
}

int clique_contains(Clique *c, int v)
{
    int k = clique_bisect(c, v);
    return (k < c->len && c->nodes[k] == v);
}

Cliques *cliques_new()
{
    Cliques *C = (Cliques *) Malloc(sizeof(Cliques), "cliques_new");
    C->len = 0;
    C->size = 0;
    C->cliques = NULL;
    C->current = NULL;
    return C;
}

Cliques *cliques_new_with_size(int size)
{
    Cliques *C = (Cliques *) Malloc(sizeof(Cliques), "cliques_new_with_size");
    C->len = 0;
    C->size = size;
    C->cliques = (Clique **) Malloc(C->size*sizeof(Clique*), "cliques_new_with_size");
    C->current = NULL;
    return C;
}

void cliques_append(Cliques *C, Clique *c)
{
    if (C->size <= C->len)
    {
        C->size = max(1, C->size * 2);
        Clique **new_cliques = (Clique **) Malloc(C->size*sizeof(Clique*), "cliques_append");
        if (C->cliques != NULL)
        {
            if (C->current) C->current = C->current - C->cliques + new_cliques;
            memcpy(new_cliques, C->cliques, C->len*sizeof(Clique*));
            free(C->cliques);
        }
        C->cliques = new_cliques;
    }
    C->cliques[C->len++] = c;
    clique_ref(c);
}

Clique *cliques_next(Cliques *C)
{
    if (!C->len)
        return NULL;

    if (!C->current)
    {
        C->current = C->cliques;
        return *C->current;
    }

    if (C->current == C->cliques + C->len - 1)
        return NULL;

    return *(++C->current);
}

void cliques_reset_iter(Cliques *C)
{
    C->current = NULL;
}

void cliques_insert(Cliques *C, Clique *c, int pos)
{
    C->cliques[pos] = c;
    clique_ref(c);
    C->len = max(C->len, pos+1);
}

void cliques_free(Cliques *C)
{
    int i;
    for (i = 0; i < C->len; ++i)
        clique_unref(C->cliques[i]);
    free(C->cliques);
    free(C);
}

#if 0
Cliques *cliques_count_sort(Cliques *C, int **count, int *len)
{
    int max_length = 0;
    int iter;
    for (iter = 0; iter < C->len; ++iter)
        max_length = max(max_length, C->cliques[iter]->len);
    int *local_count = (int *) Calloc(max_length+1, sizeof(int), "cliques_count_sort");
    for (iter = 0; iter < C->len; ++iter)
        local_count[C->cliques[iter]->len]++;
    int *first_pos = (int *) Malloc((max_length+1)*sizeof(int), "cliques_count_sort");
    first_pos[0] = 0;
    for (iter = 0; iter < max_length; ++iter)
        first_pos[iter+1] = first_pos[iter] + local_count[iter];
    Cliques *Cout = cliques_new_with_size(C->len);
    for (iter = 0; iter < C->len; ++iter)
        cliques_insert(Cout, C->cliques[iter], first_pos[C->cliques[iter]->len]++);
    *count = local_count;
    *len = max_length;
    return Cout;
}
#endif

Cliques *cliques_count_sort(Cliques *C)
{
    int max_length = 0;
    int iter;
    for (iter = 0; iter < C->len; ++iter)
        max_length = max(max_length, C->cliques[iter]->len);

    int *local_count = (int *) Calloc(max_length+1, sizeof(int), "cliques_count_sort");
    for (iter = 0; iter < C->len; ++iter)
        local_count[C->cliques[iter]->len]++;

    int *first_pos = (int *) Malloc((max_length+1)*sizeof(int), "cliques_count_sort");
    first_pos[0] = 0;
    for (iter = 0; iter < max_length; ++iter)
        first_pos[iter+1] = first_pos[iter] + local_count[iter];

    Cliques *Cout = cliques_new_with_size(C->len);
    for (iter = 0; iter < C->len; ++iter)
        cliques_insert(Cout, C->cliques[iter], first_pos[C->cliques[iter]->len]++);

    free(first_pos);
    free(local_count);

    return Cout;
}

Cliques **cliques_count_sort_array(Cliques *C, int *len)
{
    int max_length = 0;
    int iter;
    for (iter = 0; iter < C->len; ++iter)
        max_length = max(max_length, C->cliques[iter]->len);

    Cliques **Cout = (Cliques **) Malloc((max_length+1)*sizeof(Cliques *), "cliques_count_sort");
    for (iter = 0; iter < max_length+1; ++iter)
        Cout[iter] = cliques_new();

    for (iter = 0; iter < C->len; ++iter)
        cliques_append(Cout[C->cliques[iter]->len], C->cliques[iter]);

    *len = max_length;
    return Cout;
}

int cliques_contains(Cliques *C, Clique *c)
{
    int iter;
    for (iter = C->len - 1; iter > -1; --iter)
        if (C->cliques[iter] == c)
            return 1;
    return 0;
}

void get_maximal_cliques(Cliques *MC, csptr A, Clique *R, Clique *P, Clique *X)
{
    int iter;
    if (!P->len && !X->len)
    {
        cliques_append(MC, R);
        return;
    }
    Clique *potential_pivots = clique_union(P, X);
    int pivot_len = -1, pivot = -1;
    for (iter = 0; iter < potential_pivots->len; ++iter)
    {
        //~ printf("potential_pivots->nodes[iter]=%d\n", potential_pivots->nodes[iter]);
        Clique *aj = clique_new_from_nodes(A->pj[potential_pivots->nodes[iter]], A->nnzrow[potential_pivots->nodes[iter]]);
        Clique *u = clique_intersection(aj, P);
        clique_remove(u, potential_pivots->nodes[iter]);
        if (u->len > pivot_len)
        {
            pivot_len = u->len;
            pivot = potential_pivots->nodes[iter];
        }
        clique_free(aj);
        clique_free(u);
    }
    clique_free(potential_pivots);

    //~ // This can happen with an unsymmetric pattern
    //~ if (pivot == -1)
    //~ {
        //~ cliques_append(MC, R);
        //~ return;
    //~ }

    Clique *neighbours = clique_new_from_nodes(A->pj[pivot], A->nnzrow[pivot]);
    clique_remove(neighbours, pivot);
    Clique *non_neighbours = clique_exclude(P, neighbours);
    clique_free(neighbours);

    // This can happen with an unsymmetric pattern
    if (!non_neighbours->len)
    {
        cliques_append(MC, R);
        clique_free(non_neighbours);
        return;
    }

    for (iter = 0; iter < non_neighbours->len; ++iter)
    {
        clique_remove(P, non_neighbours->nodes[iter]);
        Clique *Rnew = clique_new_from_clique(R);
        clique_insert(Rnew, non_neighbours->nodes[iter]);

        Clique *Nu = clique_new_from_nodes(A->pj[non_neighbours->nodes[iter]], A->nnzrow[non_neighbours->nodes[iter]]);
        clique_remove(Nu, non_neighbours->nodes[iter]);

        Clique *Pnew = clique_intersection(P, Nu);
        Clique *Xnew = clique_intersection(X, Nu);
        clique_free(Nu);

        get_maximal_cliques(MC, A, Rnew, Pnew, Xnew);

        clique_insert(X, non_neighbours->nodes[iter]);

        clique_unref(Rnew);
        clique_unref(Pnew);
        clique_unref(Xnew);
    }
    clique_free(non_neighbours);
}

int init_blocks_cliques(csptr A, int *block_number, int **block_sizes, int **pperm, double eps)
{
    int iter, iter2;
    Cliques *MC = cliques_new();
    Clique *R = clique_new();
    Clique *P = clique_new();
    for (iter = 0; iter < A->n; ++iter)
        clique_append(P, iter);
    Clique *X = clique_new();
    get_maximal_cliques(MC, A, R, P, X);
    clique_unref(R);
    clique_unref(P);
    clique_unref(X);

    Cliques **node_to_cliques = Malloc(A->n * sizeof(Cliques*), "init_blocks_cliques");
    for (iter = 0; iter < A->n; ++iter)
        node_to_cliques[iter] = cliques_new();

    //~ for (iter = 0; iter < MC->len; ++iter)
        //~ clique_print(MC->cliques[iter]);
    printf("len mc = %d\n", MC->len);

    Clique *current;
    for (iter = 0; iter < MC->len; ++iter)
    {
        current = MC->cliques[iter];
        for (iter2 = 0; iter2 != current->len; ++iter2)
            cliques_append(node_to_cliques[current->nodes[iter2]], current);
    }

    int count_len = 0;
    Cliques **sorted_MC = cliques_count_sort_array(MC, &count_len);

    int *diag = Calloc(A->n, sizeof(int), "init_blocks_cliques");
    for (iter = 0; iter < A->n; ++iter)
        for (iter2 = 0; iter2 < A->nnzrow[iter]; ++iter2)
            if (iter == A->pj[iter][iter2])
            {
                diag[iter] = 1;
                break;
            }

    int fill, to_process_idx;
    Cliques *cliques_to_compare;
    Clique *to_compare;
    //~ Clique **current_ptr;
    //~ current_ptr = MC->cliques;
    int sorted_MC_pos = 0;
    *block_number = 0;
    while (sorted_MC_pos++ < count_len)
    {
        while ((current = cliques_next(sorted_MC[sorted_MC_pos])))
            if (!current->processed)
                break;

        if (!current)
            continue;

        if (!current->len)
        {
            current->processed = 1;
            sorted_MC_pos = 0;
            continue;
        }

        (*block_number)++;
        fill = current->len * current->len;
        for (iter = 0; iter < current->len; ++iter)
            fill -= !diag[current->nodes[iter]];
        //~ {
            //~ Clique *aj = clique_new_from_nodes(A->pj[current->nodes[iter]], A->nnzrow[current->nodes[iter]]);
            //~ fill -= !clique_contains(aj, current->nodes[iter]);
            //~ clique_free(aj);
        //~ }

        // Loop over all nodes that we have to process
        Clique *to_process = clique_new_from_clique(current);
        //~ current_to_process = to_process->nodes[0];
        to_process_idx = -1;
        while (++to_process_idx < to_process->len)
        {
            cliques_to_compare = node_to_cliques[to_process->nodes[to_process_idx]];
            // Loop over all cliques with the same node in it
            while ((to_compare = cliques_next(cliques_to_compare)))
            {
                if (to_compare == current || !to_compare->len || to_compare->processed)
                    continue;

                Clique *current_difference = clique_exclude(to_compare, current);
                int fill2 = fill + current_difference->len * to_compare->len * 2 - current_difference->len * current_difference->len;
                for (iter = 0; iter < current_difference->len; ++iter)
                    fill2 -= !diag[current_difference->nodes[iter]];
                //~ {
                    //~ Clique *aj = clique_new_from_nodes(A->pj[current_difference->nodes[iter]], A->nnzrow[current_difference->nodes[iter]]);
                    //~ fill2 -= !clique_contains(aj, current_difference->nodes[iter]);
                    //~ clique_free(aj);
                //~ }
                int total = current->len * current->len + current_difference->len * current->len * 2 + current_difference->len * current_difference->len;
                double ratio = fill2 / (double)total;
                if (ratio >= eps)
                {
                    Clique *tmp = clique_union(current, current_difference);
                    clique_copy(current, tmp);
                    clique_free(tmp);
                    clique_append_clique(to_process, current_difference);
                    fill = fill2;
                    clique_clear(to_compare);
                }
                else
                {
                    int len0 = to_compare->len;
                    clique_remove(to_compare, to_process->nodes[to_process_idx]);
                    //~ if (!len)
                        //~ reduce_count(&first_pos, count_len);
                    //~ else
                    //~ {
                    int len = to_compare->len;
                    assert(len != len0);
                    cliques_append(sorted_MC[len], to_compare);
                    //~ if (first_pos[len] != first_pos[len0] && first_pos[len0] < MC->len) //can be last in case we don't want to switch
                    //~ {
                        //~ Clique *tmp = MC->cliques[first_pos[len0]];
                        //~ MC->cliques[first_pos[len0]] = to_compare;
                        //~ MC->cliques[first_pos[len0]]
                        //~ Clique *tmp = clique_new_from_clique(MC->cliques[first_pos[len+1]]);
                        //~ clique_print(tmp);
                        //~ printf("^\n");
                        //~ clique_copy(MC->cliques[first_pos[len+1]], to_compare);
                        //~ clique_copy(to_compare, tmp);
                        //~ clique_free(tmp);
                    //~ }
                    //~ first_pos[len0]++;
                    //~ }
                }
                clique_free(current_difference);
            }
            //~ cliques_free(cliques_to_compare);
        }
        clique_free(to_process);
        sorted_MC_pos = 0;
        current->processed = 1;
    }

    //~ *block_sizes = Malloc(*block_number * sizeof(int), "init_blocks_cliques");
    int *local_block_sizes = Calloc(A->n, sizeof(int), "init_blocks_cliques");
    *pperm = Malloc(A->n * sizeof(int), "init_blocks_cliques");

    int total_len = 0;
    int *block_sizes_iter = local_block_sizes;
    int *pperm_iter = *pperm;
    int pperm_pos = 0;
    *block_number = 0;
    for (iter = 0; iter < MC->len; ++iter)
    {
        current = MC->cliques[iter];
        if (current->len)
        {
            for (iter2 = 0; iter2 < current->len; ++iter2)
                pperm_iter[current->nodes[iter2]] = iter;
        }
    }
    block_sizes_iter--;
    int prev = -1;
    for (iter = 0; iter < A->n; ++iter)
    {
        if (pperm_iter[iter] != prev)
        {
            ++(*block_number);
            ++block_sizes_iter;
        }
        (*block_sizes_iter)++;

        prev = pperm_iter[iter];
        pperm_iter[iter] = iter;
    }
    *block_sizes = Malloc(*block_number * sizeof(int), "init_blocks_cliques");
    memcpy(*block_sizes, local_block_sizes, *block_number * sizeof(int));
    free(local_block_sizes);
    //~ for (iter = 0; iter < MC->len; ++iter)
        //~ if (MC->cliques[iter]->len) clique_print(MC->cliques[iter]);

    //~ for (iter = 0; iter < MC->len; ++iter)
    //~ {
        //~ current = MC->cliques[iter];
        //~ if (current->len)
        //~ {
            //~ *block_sizes_iter = current->len;
            //~ total_len += current->len;
            //~ for (iter2 = 0; iter2 < current->len; ++iter2)
                //~ pperm_iter[current->nodes[iter2]] = pperm_pos++;
            //~ memcpy(pperm_iter,  current->nodes, current->len * sizeof(int));
            //~ pperm_iter += current->len;
            //~ block_sizes_iter++;
        //~ }
    //~ }

    printf("A->n %d, total nodes %d, nblocks %d\n", A->n, total_len, *block_number);

    cliques_free(MC);
    for (iter = 0; iter < count_len + 1; ++iter)
        cliques_free(sorted_MC[iter]);
    free(sorted_MC);
    for (iter = 0; iter < A->n; ++iter)
        cliques_free(node_to_cliques[iter]);
    free(node_to_cliques);
    free(diag);
    return 0;
}

void get_hashed_cliques(Cliques *MC, csptr A)
{
    int i, j, *ja, key, nnzrow, *ja0 = NULL;
    KeyType *group = (KeyType *)Malloc(A->n * sizeof(KeyType), "get_hashed_cliques");
    for(i = 0; i < A->n; i++)
    {
        key = 0;
        ja = A->pj[i];
        nnzrow = A->nnzrow[i];
        for(j = 0; j < nnzrow; j++)
            key += ja[j]+1;
        group[i].key = key;
        group[i].var = i;
    }
    qsort(group, A->n, sizeof(KeyType), KeyComp);

    Clique *new_clique = clique_new();
    int *iw = Calloc(A->n, sizeof(int), "get_hashed_cliques");
    int new_group = 1;
    nnzrow = 0;

    i = 0;
    while (i < A->n)
    {
        if (new_group)
        {
            // Reset iw
            for (j = 0; j < nnzrow; ++j)
                iw[ja0[j]] = 0;

            // Reset new group
            new_group = 0;

            key = group[i].key;
            nnzrow = A->nnzrow[group[i].var];

            ja0 = A->pj[group[i].var];
            for (j = 0; j < nnzrow; ++j)
                iw[ja0[j]] = 1;

            clique_unref(new_clique);
            new_clique = clique_new();
            cliques_append(MC, new_clique);
            clique_append(new_clique, group[i++].var);
            continue;
        }

        if ((new_group = (key != group[i].key)))
            continue;

        if ((new_group = (nnzrow != A->nnzrow[group[i].var])))
            continue;

        ja = A->pj[group[i].var];
        for (j = 0; j < nnzrow; ++j)
            if (iw[ja[j]] != 1)
            {
                new_group = 1;
                break;
            }

        if (new_group)
            continue;

        clique_append(new_clique, group[i++].var);
    }
    free(iw);
    clique_unref(new_clique);
}

int init_blocks_density(csptr A, int *block_number, int **block_sizes, int **pperm, double eps)
{
    int iter, iter2, iter3;
    Cliques *MC = cliques_new();
    get_hashed_cliques(MC, A);

    Cliques **node_to_cliques = Malloc(A->n * sizeof(Cliques*), "init_blocks_cliques");
    for (iter = 0; iter < A->n; ++iter)
        node_to_cliques[iter] = cliques_new();

    //~ for (iter = 0; iter < MC->len; ++iter)
        //~ clique_print(MC->cliques[iter]);
    printf("len mc = %d\n", MC->len);

    PROFILE_START("Preprocessing");

    int *ja, nnzrow;
    Clique *current;
    while ((current = cliques_next(MC)))
    {
        for (iter2 = 0; iter2 != current->len; ++iter2)
        {
            ja = A->pj[current->nodes[iter2]];
            nnzrow = A->nnzrow[current->nodes[iter2]];
            for (iter3 = 0; iter3 != nnzrow; ++iter3)
                cliques_append(node_to_cliques[ja[iter3]], current);
        }
    }
    cliques_reset_iter(MC);

    Cliques *sorted_MC = cliques_count_sort(MC);

    int fill, to_process_idx;
    Cliques *cliques_to_compare;
    Clique *to_compare;
    *block_number = 0;
    //~ Cliques totals = cliques_new();
    for (iter = 0; iter < MC->len; ++iter)
    {
        current = MC->cliques[iter];
        Clique *total = clique_new();
        for (iter2 = 0; iter2 < current->len; ++iter2)
        {
            Clique *tmp = clique_new_from_nodes(A->pj[current->nodes[iter2]], A->nnzrow[current->nodes[iter2]]);
            Clique *tmp2 = clique_new_from_clique(total);
            clique_free(total);
            total = clique_union(tmp, tmp2);
            clique_free(tmp);
            clique_free(tmp2);
        }
        clique_add_neighbours(current, total);
        Clique *tmp = clique_new_from_clique(total);
        clique_add_neighbour_blocks(current, tmp);
        clique_unref(total);
        clique_unref(tmp);
    }

    int *node_to_block_size = Calloc(A->n, sizeof(int), "init_blocks_cliques");
    for (iter = 0; iter < MC->len; ++iter)
    {
        current = MC->cliques[iter];
        for (iter2 = 0; iter2 < current->len; ++iter2)
            node_to_block_size[current->nodes[iter2]] = current->len;
    }

    Clique **node_to_clique = Malloc(A->n * sizeof(Clique*), "init_blocks_cliques");
    for (iter = 0; iter < MC->len; ++iter)
    {
        current = MC->cliques[iter];
        for (iter2 = 0; iter2 < current->len; ++iter2)
            node_to_clique[current->nodes[iter2]] = current;
    }
    PROFILE_END("Preprocessing");

    //~ while ((current = cliques_next(sorted_MC)))
    for (iter = 0; iter < MC->len; ++iter)
    {
        current = MC->cliques[iter];
        if (!current || current->processed || !current->len)
            continue;

        (*block_number)++;

        // Compute the amount of nonzeros in the clique
        fill = 0;
        for (iter2 = 0; iter2 < current->len; ++iter2)
            fill += A->nnzrow[current->nodes[iter2]];
        assert(fill == ((Clique *)current->neighbours)->len * current->len);
        // Compensate for symmetry
        fill = 2 * fill - current->len * current->len;
//~ 
        //~ Clique *neighbour_blocks = clique_new();
        //~ Clique *prev_ptr = NULL;
        //~ for (iter2 = 0; iter2 < ((Clique *)current->neighbours)->len; ++iter2)
        //~ {
            //~ if (prev_ptr == node_to_clique[((Clique *)current->neighbours)->nodes[iter2]])
                //~ continue;
//~ 
            //~ prev_ptr = node_to_clique[((Clique *)current->neighbours)->nodes[iter2]];
            //~ if (clique_contains(neighbour_blocks, ((Clique *)current->neighbours)->nodes[iter2]))
                //~ continue;
//~ 
            //~ Clique *tmp = clique_union(neighbour_blocks, node_to_clique[((Clique *)current->neighbours)->nodes[iter2]]);
            //~ clique_free(neighbour_blocks);
            //~ neighbour_blocks = tmp;
        //~ }

        // Loop over all nodes that we have to process
        Clique *to_process = clique_new_from_clique(current);
        Cliques *cliques_to_compare_done = cliques_new();
        to_process_idx = -1;
        while (++to_process_idx < to_process->len)
        {
            cliques_to_compare = node_to_cliques[to_process->nodes[to_process_idx]];
            // Loop over all cliques with the same node in it
            cliques_reset_iter(cliques_to_compare);
            while ((to_compare = cliques_next(cliques_to_compare)))
            {
                if (to_compare == current || !to_compare->len || to_compare->processed || cliques_contains(cliques_to_compare_done, to_compare))
                    continue;

                Clique *current_neighbour_blocks = clique_union((Clique *)current->neighbour_blocks, (Clique *)to_compare->neighbour_blocks);
                Clique *current_difference = clique_exclude((Clique *)to_compare->neighbours, current);
                int fill2 = fill + 2 * to_compare->len * current_difference->len - to_compare->len * to_compare->len;

                int total_len = current->len + to_compare->len;
                //~ Clique *total_neighbours = clique_new_from_clique(current_neighbours);
                //~ Clique *total_blocks = clique_new_from_clique(neighbour_blocks);
                //~ prev_ptr = NULL;
                //~ for (iter2 = 0; iter2 < ((Clique *)to_compare->neighbours)->len; ++iter2)
                //~ {
                    //~ if (prev_ptr == node_to_clique[((Clique *)to_compare->neighbours)->nodes[iter2]])
                        //~ continue;
//~ 
                    //~ prev_ptr = node_to_clique[((Clique *)to_compare->neighbours)->nodes[iter2]];
                    //~ if (clique_contains(total_blocks, ((Clique *)to_compare->neighbours)->nodes[iter2]))
                        //~ continue;
//~ 
                    //~ Clique *tmp = clique_union(total_blocks, node_to_clique[((Clique *)to_compare->neighbours)->nodes[iter2]]);
                    //~ clique_free(total_blocks);
                    //~ total_blocks = tmp;
                //~ }
                //~ int total_nodes = 2 * total_blocks->len * total_len - total_len * total_len;
                //~ int total_nodes = 2 * current_neighbours->len * total_len - total_len * total_len;
                int total_nodes = 2 * current_neighbour_blocks->len * total_len - total_len * total_len;
                double ratio = fill2 / (double) total_nodes;
                //~ printf("fill2: %d, total: %d, ratio: %f\n", fill2, total_nodes, ratio);
                //~ printf("fill: %d, total: %d, ratio: %f\n", fill, 2 * total->len * current->len - current->len * current->len, fill / (double) ( 2 * total->len * current->len - current->len * current->len));
                if (ratio >= eps)
                {
                    Clique *current_neighbours = clique_union((Clique *)current->neighbours, (Clique *)to_compare->neighbours);
                    Clique *tmp = clique_union(current, to_compare);
                    assert(tmp->len == total_len);
                    clique_copy(current, tmp);
                    clique_free(tmp);
                    clique_append_clique(to_process, to_compare);
                    fill = fill2;
                    clique_add_neighbours(current, current_neighbours);
                    clique_add_neighbour_blocks(current, current_neighbour_blocks);

                    PROFILE_START("Blocks Update");
                    Cliques *had = cliques_new();
                    for (iter = 0; iter < current->len; ++iter)
                    {
                        Cliques *cliques_to_update = node_to_cliques[current->nodes[iter]];
                        for (iter2 = 0; iter2 < cliques_to_update->len; ++iter2)
                        {
                            if (cliques_to_update->cliques[iter2]->processed || cliques_contains(had, cliques_to_update->cliques[iter2]))
                                continue;

                            cliques_append(had, cliques_to_update->cliques[iter2]);
                            Clique *tmp = clique_union(cliques_to_update->cliques[iter2]->neighbour_blocks, current);
                            clique_add_neighbour_blocks(cliques_to_update->cliques[iter2], tmp);
                            clique_unref(tmp);
                        }
                    }
                    cliques_free(had);
                    PROFILE_END("Blocks Update");
                    clique_unref(current_neighbours);
//~ 
                    //~ for (iter2 = 0; iter2 < ((Clique *)to_compare->neighbours)->len; ++iter2)
                        //~ node_to_clique[((Clique *)to_compare->neighbours)->nodes[iter2]] = current;
                    //~ clique_free(neighbour_blocks);
                    //~ neighbour_blocks = total_blocks;
                    //~ for (iter2 = 0; iter2 < to_compare->len; ++iter2)
                        //~ node_to_clique[to_compare->nodes[iter2]] = current;

                    clique_clear(to_compare);
                    to_compare->processed = 1;
                }
                //~ else
                //~ {
                    //~ clique_free(total_blocks);
                //~ }
                clique_unref(current_neighbour_blocks);
                clique_free(current_difference);
                cliques_append(cliques_to_compare_done, to_compare);
            }
        }
        cliques_free(cliques_to_compare_done);
        clique_free(to_process);
        current->processed = 1;
    }

    FILE *fp;
    fp = fopen("profile1.txt", "w");
    SAVE_PROFILES(fp);
    fclose(fp);

    int *local_block_sizes = Calloc(A->n, sizeof(int), "init_blocks_cliques");
    *pperm = Malloc(A->n * sizeof(int), "init_blocks_cliques");

    int total_len = 0;
    int *block_sizes_iter = local_block_sizes;
    int *pperm_iter = *pperm;
    int pperm_pos = 0;
    for (iter = 0; iter < MC->len; ++iter)
    {
        current = MC->cliques[iter];
        if (current->len)
        {
            *block_sizes_iter = current->len;
            total_len += current->len;
            for (iter2 = 0; iter2 < current->len; ++iter2)
                pperm_iter[current->nodes[iter2]] = pperm_pos++;
            block_sizes_iter++;
        }
    }

    *block_sizes = Malloc(*block_number * sizeof(int), "init_blocks_cliques");
    memcpy(*block_sizes, local_block_sizes, *block_number * sizeof(int));
    free(local_block_sizes);

    printf("A->n %d, total nodes %d, nblocks %d MC->len %d\n", A->n, total_len, *block_number, MC->len);

    cliques_free(MC);
    cliques_free(sorted_MC);
    for (iter = 0; iter < A->n; ++iter)
        cliques_free(node_to_cliques[iter]);
    free(node_to_cliques);
    return 0;
}
