#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "globheads.h"
#include "protos.h"

//~ #define PROFILE_ON
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

typedef struct clique{
    int len; // amount of nodes
    int *nodes;
    int size; //size of nodes
    int *end; //end ptr
    int ref; //ref count
    int processed;
    int fill;
    void *neighbours;
    void *neighbour_blocks;
} Clique;

typedef struct cliques{
    int len; // amount of cliques
    int size; // size of cliques
    Clique **cliques;
    Clique **end;
} Cliques;

static Clique *clique_new()
{
    PROFILE_START("New");
    Clique *c = (Clique *) Malloc(sizeof(Clique), "clique_new");
    c->len = 0;
    c->size = 0;
    c->nodes = NULL;
    c->end = NULL;
    c->ref = 1;
    c->processed = 0;
    c->fill = 0;
    c->neighbours = NULL;
    c->neighbour_blocks = NULL;
    PROFILE_END("New");
    return c;
}

static Clique *clique_new_with_size(int size)
{
    Clique *c = (Clique *) Malloc(sizeof(Clique), "clique_new_with_size");
    c->len = 0;
    c->size = size;
    c->nodes = (int *) Malloc(c->size*sizeof(int), "clique_new_with_size");
    c->end = c->nodes;
    c->ref = 1;
    c->processed = 0;
    c->fill = 0;
    c->neighbours = NULL;
    c->neighbour_blocks = NULL;
    return c;
}

static Clique *clique_new_from_nodes_ref(int *nodes, int len)
{
    PROFILE_START("New from nodes ref");
    Clique *c = (Clique *) Malloc(sizeof(Clique), "clique_new_from_nodes_ref");
    c->len = len;
    c->size = len;
    c->nodes = nodes;
    c->end = c->nodes + len;
    c->ref = 1;
    c->processed = 0;
    c->fill = 0;
    c->neighbours = NULL;
    c->neighbour_blocks = NULL;
    PROFILE_END("New from nodes ref");
    return c;
}

static void clique_unref(Clique *c);
static void clique_remove_neighbours(Clique * const c);

static void clique_unset(Clique *c)
{
    if (c->nodes)
        free(c->nodes);
    if (c->neighbours)
        clique_remove_neighbours(c);
    if (c->neighbour_blocks)
        clique_unref((Clique *) c->neighbour_blocks);
}

static void clique_free(Clique *c)
{
    clique_unset(c);
    free(c);
}

static void clique_ref(Clique * const c)
{
    c->ref++;
}

static void clique_unref(Clique *c)
{
    if (!--c->ref)
        clique_free(c);
}

static void clique_remove_neighbours(Clique * const c)
{
    if (!--((Clique *)c->neighbours)->ref)
        free((Clique *)c->neighbours);
    c->neighbours = NULL;
}

static void clique_add_neighbours(Clique * const c, Clique * const n)
{
    if (c->neighbours)
        clique_remove_neighbours(c);

    c->neighbours = (void *)n;
    clique_ref(n);
}

static void clique_remove_neighbour_blocks(Clique * const c)
{
    clique_unref((Clique *)c->neighbour_blocks);
    c->neighbour_blocks = NULL;
}

static void clique_add_neighbour_blocks(Clique * const c, Clique * const n)
{
    if (c->neighbour_blocks)
        clique_unref((Clique *)c->neighbour_blocks);

    c->neighbour_blocks = (void *)n;
    clique_ref(n);
}

static void clique_set_size(Clique * const c, int size)
{
    if (c->size < size)
    {
        free(c->nodes);
        c->size = size;
        c->nodes = (int *) Malloc(c->size*sizeof(int), "clique_set_size");
    }
}

static void clique_copy(Clique * const dest, Clique * const src)
{
    clique_set_size(dest, src->size);

    if (dest->neighbours)
        clique_remove_neighbours(dest);

    if (dest->neighbour_blocks)
        clique_remove_neighbour_blocks(dest);

    dest->len = src->len;
    dest->size = src->size;
    dest->processed = src->processed;
    dest->fill = src->fill;
    memcpy(dest->nodes, src->nodes, dest->len*sizeof(int));
    dest->end = dest->nodes + dest->len;

    if (src->neighbours)
        clique_add_neighbours(dest, (Clique *)src->neighbours);

    if (src->neighbour_blocks)
        clique_add_neighbour_blocks(dest, (Clique *)src->neighbour_blocks);
}

static Clique *clique_new_from_clique(Clique * const src)
{
    Clique *dest = clique_new();
    clique_copy(dest, src);
    return dest;
}

void clique_clear(Clique * const c)
{
    if (c->nodes)
        free(c->nodes);

    c->nodes = NULL;
    c->len = 0;
    c->size = 0;
    c->end = NULL;
}

static void clique_move(Clique *dest, Clique *src)
{
    memcpy(dest, src, sizeof(Clique));
    src->nodes = NULL;
    src->neighbours = NULL;
    src->neighbour_blocks = NULL;
}

static void clique_set_from_clique(Clique * const dest, Clique * const src)
{
    if (dest->size < src->size)
    {
        clique_copy(dest, src);
        return;
    }

    dest->len = src->len;
    memcpy(dest->nodes, src->nodes, src->len*sizeof(int));
    dest->end = dest->nodes + dest->len;
}

static void clique_append(Clique * const c, int v)
{
    if (c->size <= c->len)
    {
        PROFILE_START("Realloc clique");
        c->size = max(MAX_BLOCK_SIZE, c->size * 2);
        int *new_clique = (int *) Malloc(c->size*sizeof(int), "clique_append");
        if (c->nodes != NULL)
        {
            memcpy(new_clique, c->nodes, c->len*sizeof(int));
            free(c->nodes);
        }
        c->nodes = new_clique;
        c->end = c->nodes + c->len;
        PROFILE_END("Realloc clique");
    }
    c->len++;
    *c->end++ = v;
}

static void clique_append_unsafe(Clique * const c, int v)
{
    *c->end++ = v;
}

static void clique_append_clique(Clique * const x, Clique * const y)
{
    if (!y->len)
        return;

    if (x->size < x->len + y->len)
    {
        x->size = x->len + y->len;
        int *new_clique = (int *) Malloc(x->size*sizeof(int), "clique_append_clique");
        if (x->nodes != NULL)
        {
            memcpy(new_clique, x->nodes, x->len*sizeof(int));
            free(x->nodes);
        }
        x->nodes = new_clique;
        x->end = x->nodes + x->len;
    }
    memcpy(x->end, y->nodes, y->len*sizeof(int));
    x->len += y->len;
    x->end += y->len;
}

void clique_print(Clique * const c)
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

static Clique *clique_union(Clique * const x, Clique * const y)
{
    PROFILE_START("Union");

    Clique *out = clique_new_with_size(x->len + y->len);

    int *a = x->nodes;
    int *b = y->nodes;

    int last = -1;
    int c;

    while (a != x->end && b != y->end)
    {
        c = *a < *b ? *a++ : *b++;
        if (c != last)
            clique_append_unsafe(out, last = c);
    }

    if (a != x->end && *a == last) a++;
    if (a != x->end)
    {
        const int len = x->end - a;
        memcpy(out->end, a, sizeof(int)*len);
        out->end += len;
    }
    if (b != y->end && *b == last) b++;
    if (b != y->end)
    {
        const int len = y->end - b;
        memcpy(out->end, b, sizeof(int)*len);
        out->end += len;
    }
    out->len = out->end - out->nodes;
    PROFILE_END("Union");
    return out;
}

static void clique_set_union(Clique * const out, Clique * const x, Clique * const y)
{
    PROFILE_START("Set Union");

    clique_set_size(out, x->len + y->len);

    int *a = x->nodes;
    int *b = y->nodes;

    int last = -1;
    int c;
    out->end = out->nodes;

    while (a != x->end && b != y->end)
    {
        c = *a < *b ? *a++ : *b++;
        if (c != last)
            clique_append_unsafe(out, last = c);
    }

    if (a != x->end && *a == last) a++;
    if (a != x->end)
    {
        const int len = x->end - a;
        memcpy(out->end, a, sizeof(int)*len);
        out->end += len;
    }
    if (b != y->end && *b == last) b++;
    if (b != y->end)
    {
        const int len = y->end - b;
        memcpy(out->end, b, sizeof(int)*len);
        out->end += len;
    }
    out->len = out->end - out->nodes;
    PROFILE_END("Set Union");
}

static int clique_exclude_len(Clique * const x, Clique * const y)
{
    PROFILE_START("Exclude len");
    int len = x->len;

    int *a = x->nodes;
    int *b = y->nodes;

    while (a != x->end && b != y->end)
    {
        if (*a < *b) a++;
        else if (*a == *b++)
            len--;
    }

    PROFILE_END("Exclude len");
    return len;
}

static int clique_union_len(Clique * const x, Clique * const y)
{
    return clique_exclude_len(x, y) + y->len;
}

static Cliques *cliques_new()
{
    PROFILE_START("Cliques New");
    Cliques *C = (Cliques *) Malloc(sizeof(Cliques), "cliques_new");
    C->len = 0;
    C->size = 0;
    C->cliques = NULL;
    C->end = NULL;
    PROFILE_END("Cliques New");
    return C;
}

static void cliques_set_size(Cliques *C, int size)
{
    if (C->size < size)
    {
        C->size = max(size, C->size * 2);
        Clique **new_cliques = (Clique **) Malloc(C->size*sizeof(Clique*), "cliques_set_size");
        if (C->cliques != NULL)
        {
            memcpy(new_cliques, C->cliques, C->len*sizeof(Clique*));
            free(C->cliques);
        }
        C->cliques = new_cliques;
        C->end = C->cliques + C->len;
    }
}

static void cliques_free_nodes(Cliques *C)
{
    int i;
    for (i = 0; i < C->len; ++i)
        clique_unset(C->cliques[i]);
    free(C->cliques);
    free(C);
}

static void cliques_append_noref(Cliques *C, Clique *c)
{
    if (C->size <= C->len)
    {
        C->size = max(MAX_BLOCK_SIZE, C->size * 2);
        Clique **new_cliques = (Clique **) Malloc(C->size*sizeof(Clique*), "cliques_append");
        if (C->cliques != NULL)
        {
            memcpy(new_cliques, C->cliques, C->len*sizeof(Clique*));
            free(C->cliques);
        }
        C->cliques = new_cliques;
        C->end = new_cliques + C->len;
    }
    C->len++;
    *C->end++ = c;
}

static void cliques_append_noref_unsafe(Cliques *C, Clique *c)
{
    *C->end++ = c;
}

static void cliques_remove_noref(Cliques *C, Clique *c)
{
    PROFILE_START("Cliques Remove");
    int found = 0;

    // Do one bisection step. More hurts performance because the next loop
    // is very fast
    Clique **iter = C->cliques + C->len / 2;
    if (c < *iter)
        iter = C->cliques;

    while (iter != C->end)
        if (c == *iter++)
        {
            found = 1;
            break;
        }
    if (found)
    {
        memmove(iter-1, iter, (C->end-- -iter) * sizeof(Clique *));
        C->len--;
    }

    PROFILE_END("Cliques Remove");
}

static void cliques_union_noref(Cliques *out, Cliques *X, Cliques *Y)
{
    PROFILE_START("Cliques Union Noref");

    Clique **a = X->cliques;
    Clique **b = Y->cliques;
    Clique *merge;
    Clique *last = NULL;

    cliques_set_size(out, X->len + Y->len);

    out->end = out->cliques;

    if (!X->len)
    {
        while (b != Y->end)
        {
            merge = *b++;
            if (!merge->processed && last != merge)
                cliques_append_noref_unsafe(out, last = merge);
        }
        out->len = out->end - out->cliques;
        PROFILE_END("Cliques Union Noref");
        return;
    }

    while (a != X->end && b != Y->end)
    {
        merge = *a < *b ? *a++ : *b++;
        if (!merge->processed && last != merge)
            cliques_append_noref_unsafe(out, last = merge);
    }

    while (a != X->end)
    {
        merge = *a++;
        if (!merge->processed && last != merge)
            cliques_append_noref_unsafe(out, last = merge);
    }

    while (b != Y->end)
    {
        merge = *b++;
        if (!merge->processed && last != merge)
            cliques_append_noref_unsafe(out, last = merge);
    }

    out->len = out->end - out->cliques;

    PROFILE_END("Cliques Union Noref");
}

static void cliques_union_noref_inc_processed(Cliques *out, Cliques *X, Cliques *Y)
{
    PROFILE_START("Cliques Union Noref Processed");

    Clique **a = X->cliques;
    Clique **b = Y->cliques;
    Clique *merge;
    Clique *last = NULL;

    cliques_set_size(out, X->len + Y->len);

    out->end = out->cliques;

    if (!X->len)
    {
        while (b != Y->end)
        {
            merge = *b++;
            if (merge->len && last != merge)
                cliques_append_noref_unsafe(out, last = merge);
        }
        out->len = out->end - out->cliques;
        PROFILE_END("Cliques Union Noref Processed");
        return;
    }

    while (a != X->end && b != Y->end)
    {
        merge = *a < *b ? *a++ : *b++;
        if (merge->len && merge != last)
            cliques_append_noref_unsafe(out, last = merge);
    }

    while (a != X->end)
    {
        merge = *a++;
        if (merge->len && merge != last)
            cliques_append_noref_unsafe(out, last = merge);
    }

    while (b != Y->end)
    {
        merge = *b++;
        if (merge->len && merge != last)
            cliques_append_noref_unsafe(out, last = merge);
    }

    out->len = out->end - out->cliques;

    //~ a = out->cliques - 1;
    //~ while ((++a) != out->end)
    //~ printf("%p ", *a);
    //~ printf("\n");

    PROFILE_END("Cliques Union Noref Processed");
}

static void cliques_clear_noref(Cliques *C)
{
    C->len = 0;
    C->end = C->cliques;
}

static void cliques_array_free_noref(Cliques *C)
{
    free(C->cliques);
}

static void cliques_free_noref(Cliques *C)
{
    free(C->cliques);
    free(C);
}

static int cliques_first_sort_fun(const void *vx, const void *vy)
{
    Clique *x = *(Clique **)vx;
    Clique *y = *(Clique **)vy;
    if (!x->len)
        return -1;
    if (!y->len)
        return 1;
    if (x->nodes[0] < y->nodes[0])
        return -1;
    return 1;
}

static void cliques_first_sort(Cliques *C)
{
    qsort(C->cliques, C->len, sizeof(Clique *), cliques_first_sort_fun);
}

static void get_hashed_cliques(Cliques *MC, csptr A)
{
    PROFILE_START("Saad part");
    int i, i2, j, *ja, key, nnzrow, *ja0 = NULL;
    KeyType *group = (KeyType *) Malloc(A->n * sizeof(KeyType), "get_hashed_cliques");
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
    PROFILE_END("Saad part");

    PROFILE_START("My part");
    Clique *new_clique;
    int *iw = (int *) Calloc(A->n, sizeof(int), "get_hashed_cliques");
    int *processed = (int *) Calloc(A->n, sizeof(int), "get_hashed_cliques");
    nnzrow = 0;

    Clique **clique_iter = MC->cliques;
    for (i = 0; i < A->n; ++i)
    {
        if (processed[i])
            continue;

        // Reset iw
        for (j = 0; j < nnzrow; ++j)
            iw[ja0[j]] = 0;

        key = group[i].key;
        nnzrow = A->nnzrow[group[i].var];

        ja0 = A->pj[group[i].var];
        for (j = 0; j < nnzrow; ++j)
            iw[ja0[j]] = 1;

        new_clique = *clique_iter++;
        clique_append(new_clique, group[i].var);

        for (i2 = i + 1; i2 < A->n; ++i2)
        {
            if (processed[i2])
                continue;

            if (key != group[i2].key)
                break;

            if (nnzrow != A->nnzrow[group[i2].var])
                continue;

            int new_group = 0;
            j = nnzrow;
            ja = A->pj[group[i2].var];
            while (j && (new_group = (!iw[ja[--j]])));

            if (new_group)
                continue;

            clique_append(new_clique, group[i2].var);
            processed[i2] = 1;
        }
    }
    MC->len = clique_iter - MC->cliques;
    MC->end = clique_iter;
    free(iw);
    free(group);
    free(processed);
    PROFILE_END("My part");
}

int init_blocks_density(csptr A, int *block_number, int **block_sizes, int **pperm, double eps)
{
    PROFILE_START("Hashed cliques");
    Clique **MC_iter;
    int iter, iter2;

    Cliques *tmp_MC = cliques_new();
    Clique *tmp_MC_cliques = (Clique *) Calloc(A->n, sizeof(Clique), "init_blocks_density");
    Clique *tmp_MC_cliques_iter = tmp_MC_cliques;
    cliques_set_size(tmp_MC, A->n);
    for (iter = 0; iter < A->n; iter++)
        cliques_append_noref_unsafe(tmp_MC, tmp_MC_cliques_iter++);

    get_hashed_cliques(tmp_MC, A);
    PROFILE_END("Hashed cliques");

    PROFILE_START("Resort");
    // To keep as much as possible the structure of the original matrix
    cliques_first_sort(tmp_MC);
    PROFILE_END("Resort");

    PROFILE_START("MC copy");
    Cliques *MC = cliques_new();
    Clique *MC_cliques = (Clique *) Calloc(tmp_MC->len, sizeof(Clique), "init_blocks_density");
    Clique *MC_cliques_iter = MC_cliques;
    MC_iter = tmp_MC->cliques;
    cliques_set_size(MC, tmp_MC->len);
    for (MC_iter = tmp_MC->cliques; MC_iter != tmp_MC->end; ++MC_iter)
    {
        clique_move(MC_cliques_iter, *MC_iter);
        cliques_append_noref_unsafe(MC, MC_cliques_iter++);
    }
    MC->len = MC->end - MC->cliques;

    cliques_free_nodes(tmp_MC);
    free(tmp_MC_cliques);

    PROFILE_END("MC copy");

    PROFILE_START("Preprocessing 1");

    //Fast allocation of all cliques
    Cliques **node_to_cliques = Malloc(A->n * sizeof(Cliques*), "init_blocks_density");
    Cliques *cliques_array = Calloc(A->n, sizeof(Cliques), "init_blocks_density");
    for (iter = 0; iter < A->n; ++iter)
        node_to_cliques[iter] = cliques_array++;

    printf("len mc = %d\n", MC->len);

    PROFILE_END("Preprocessing 1");

    PROFILE_START("Preprocessing 2");
    int *ja, nnzrow;
    Clique *current;
    for (MC_iter = MC->cliques; MC_iter != MC->end; ++MC_iter)
    {
        current = *MC_iter;
        // Make a nodes > cliques map
        ja = A->pj[*current->nodes];
        nnzrow = A->nnzrow[*current->nodes];
        PROFILE_START("Append part");
        for (iter2 = 0; iter2 != nnzrow; ++iter2)
            //~ cliques_append_noref(node_to_cliques[ja[iter2]], current);
            cliques_append_noref(node_to_cliques[ja[iter2]], current);
        PROFILE_END("Append part");

        // Set the neighbours of every clique
        Clique * const total = clique_new_from_nodes_ref(ja, nnzrow);
        clique_add_neighbours(current, total);

        // Need a new clique because neighbours and neighbour_blocks will be different later
        Clique * const total2 = clique_new_from_clique(total);
        clique_add_neighbour_blocks(current, total2);

        current->fill = (2 * total->len - current->len) * current->len;

        clique_unref(total);
        clique_unref(total2);
    }
    PROFILE_END("Preprocessing 2");

    int *clique_iter, to_process_iter;
    Clique *to_compare;
    *block_number = 0;

    Clique * const to_process = clique_new();
    Clique * const new_current = clique_new();

    Cliques *cliques_tmp = cliques_new();
    Cliques *cliques_tmp_ptr;
    Cliques *cliques_to_compare = cliques_new();
    Cliques *cliques_to_update = cliques_new();
    Cliques * const new_neighbour_blocks_list = cliques_new();

    Clique **clique_to_update_iter;
    Clique **clique_to_update_iter2;

    PROFILE_START("Main loop");
    for (MC_iter = MC->cliques; MC_iter != MC->end; ++MC_iter)
    {
        current = *MC_iter;
        if (!current || current->processed || !current->len)
            continue;

        (*block_number)++;

        clique_set_from_clique(to_process, current);
        to_process_iter = 0;

        // Start comparing at 1 and set the first element to current so we can skip that one and don't have to check for it
        int to_compare_iter = 0;
        //~ cliques_append_noref(cliques_to_compare, current);

        // Loop over all nodes that we have to process
        while (to_process_iter + to_process->nodes != to_process->end)
        {
            to_compare_iter--;
            // Get the unique cliques that we need to compare to. The merge function excludes processed cliques
            PROFILE_START("Merge to compare");
            while (to_process_iter + to_process->nodes != to_process->end)
            {
                cliques_union_noref(cliques_tmp, cliques_to_compare, node_to_cliques[*(to_process->nodes + to_process_iter++)]);
                cliques_tmp_ptr = cliques_tmp;
                cliques_tmp = cliques_to_compare;
                cliques_to_compare = cliques_tmp_ptr;
            }
            cliques_remove_noref(cliques_to_compare, current);
            PROFILE_END("Merge to compare");
            // Loop over all cliques with the same node in it
            PROFILE_START("Compare");
            while (++to_compare_iter < cliques_to_compare->len)
            {
                to_compare = cliques_to_compare->cliques[to_compare_iter];

                PROFILE_START("Real compare");
                const int union_len = clique_union_len((Clique *)current->neighbour_blocks, (Clique *)to_compare->neighbour_blocks);

                // Compute the off-diagonal elements that will be added to be able to compute the real fill
                const int exclude_len = clique_exclude_len((Clique *)to_compare->neighbours, current);
                const int fill2 = current->fill + (2 * exclude_len - to_compare->len) * to_compare->len;

                const int total_len = current->len + to_compare->len;
                const int total_nodes = (2 * union_len - total_len) * total_len;
                const double ratio = fill2 / (double)total_nodes;
                if (ratio < eps)
                {
                    PROFILE_END("Real compare");
                    continue;
                }

                PROFILE_START("Yay we need to merge");

                clique_set_union(new_current, current, to_compare);

                PROFILE_START("Blocks Test");

                for (clique_iter = new_current->nodes; clique_iter != new_current->end; ++clique_iter)
                {
                    cliques_union_noref_inc_processed(cliques_tmp, cliques_to_update, node_to_cliques[*clique_iter]);
                    cliques_tmp_ptr = cliques_tmp;
                    cliques_tmp = cliques_to_update;
                    cliques_to_update = cliques_tmp_ptr;
                }
                cliques_remove_noref(cliques_to_update, current);
                cliques_remove_noref(cliques_to_update, to_compare);

                Clique *clique_to_update;
                double ratio2 = 1.0;

                while (cliques_to_update->len > new_neighbour_blocks_list->len)
                    cliques_append_noref(new_neighbour_blocks_list, clique_new());
                clique_to_update_iter2 = new_neighbour_blocks_list->cliques;

                // Check whether we don't influence other cliques too much by this operation
                for (clique_to_update_iter = cliques_to_update->cliques; clique_to_update_iter != cliques_to_update->end; ++clique_to_update_iter)
                {
                    clique_to_update = *clique_to_update_iter;
                    Clique *tmp = *clique_to_update_iter2++;
                    clique_set_union(tmp, (Clique *)clique_to_update->neighbour_blocks, new_current);
                    ratio2 = clique_to_update->fill / (double)(( 2 * tmp->len - clique_to_update->len) * clique_to_update->len);
                    if (ratio2 < eps)
                        break;
                }

                // And stop if we do
                if (ratio2 < eps)
                {
                    cliques_clear_noref(cliques_to_update);
                    PROFILE_END("Blocks Test");
                    PROFILE_END("Yay we need to merge");
                    PROFILE_END("Real compare");
                    continue;
                }
                PROFILE_END("Blocks Test");
                PROFILE_START("Blocks Update");

                clique_to_update_iter2 = new_neighbour_blocks_list->cliques;
                for (clique_to_update_iter = cliques_to_update->cliques; clique_to_update_iter != cliques_to_update->end; ++clique_to_update_iter)
                    clique_set_from_clique((Clique *)((*clique_to_update_iter)->neighbour_blocks), *clique_to_update_iter2++);

                cliques_clear_noref(cliques_to_update);
                PROFILE_END("Blocks Update");

                Clique * const new_neighbour_blocks = clique_union((Clique *)current->neighbour_blocks, (Clique *)to_compare->neighbour_blocks);

                assert(new_current->len == total_len);
                clique_set_from_clique(current, new_current);

                clique_append_clique(to_process, to_compare);
                current->fill = fill2;

                clique_add_neighbour_blocks(current, new_neighbour_blocks);

                clique_unref(new_neighbour_blocks);

                to_compare->len = 0;
                to_compare->processed = 1;
                PROFILE_END("Yay we need to merge");
                PROFILE_END("Real compare");
            }
            PROFILE_END("Compare");
        }
        cliques_clear_noref(cliques_to_compare);
        current->processed = 1;
    }

    clique_free(to_process);
    clique_free(new_current);

    cliques_free_noref(cliques_to_compare);
    cliques_free_noref(cliques_to_update);
    cliques_free_noref(new_neighbour_blocks_list);

    PROFILE_END("Main loop");

    PROFILE_START("Postprocessing");

    *block_sizes = (int *) Malloc(*block_number * sizeof(int), "init_blocks_density");
    *pperm = (int *) Malloc(A->n * sizeof(int), "init_blocks_density");

    //~ *block_sizes = (int *) Malloc(A->n * sizeof(int), "init_blocks_density");
    //~ for (iter=0;iter<A->n;++iter)
    //~ {
        //~ (*pperm)[iter] = iter;
        //~ (*block_sizes)[iter] = 1;
        //~ (*block_number) = A->n;
    //~ }

    int total_len = 0;
    int *block_sizes_iter = *block_sizes;
    int *pperm_iter = *pperm;
    int pperm_pos = 0;
    for (MC_iter = MC->cliques; MC_iter != MC->end; ++MC_iter)
    {
        current = *MC_iter;
        if (current->len)
        {
            *block_sizes_iter = current->len;
            total_len += current->len;
            for (clique_iter = current->nodes; clique_iter != current->end; ++clique_iter)
                pperm_iter[*clique_iter] = pperm_pos++;
            block_sizes_iter++;
        }
    }

    printf("A->n %d, total nodes %d, nblocks %d nblocks %ld MC->len %d\n", A->n, total_len, *block_number, block_sizes_iter - *block_sizes, MC->len);

    Cliques **cliques_iter = node_to_cliques;
    Cliques **cliques_end = node_to_cliques + A->n;
    for (cliques_iter = node_to_cliques; cliques_iter != cliques_end; ++cliques_iter)
        cliques_array_free_noref(*cliques_iter);
    free(*node_to_cliques);
    free(node_to_cliques);
    cliques_free_nodes(MC);
    free(MC_cliques);
    PROFILE_END("Postprocessing");

    SAVE_PROFILES("profile1.txt");

    return 0;
}
