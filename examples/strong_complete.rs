// Code to prove that if the square of the linegraph of a graph G with maximum degree Δ
// is a clique, then G has at most 5/4 Δ^2 edges

extern crate flag_algebra;
extern crate local_flags;

use flag_algebra::*;
use flag_algebra::flags::{Graph, Connected};
use local_flags::Degree;

type F = SubClass<Graph, Connected>;
type N = i64;
type V = QFlag<N, F>;

fn strong_degree() -> V {
    let edge: F = Graph::new(2, &[(0, 1)]).into();
    let t = Type::from_flag(&edge);
    Basis::new(4).with_type(t).from_indicator(
        | g: &F, _ | g.content.edge(2,3)                                         
    )
}

fn one(n: usize) -> V {
    let single_vertex: F = Graph::new(1, &[]).into();
    single_vertex.project(n).untype()
}

// Sum of graphs containing a 2K2
fn not_strong_clique(size: usize) -> V {
    Basis::new(size).from_indicator( |g: &F, _|{
        let n = g.content.size();
        for (a1, a2) in g.content.edges() {
            for b1 in 0..n {
                if b1 != a1 && b1 != a2 {
                    for b2 in (b1 + 1)..n {
                        if g.content.edge(b1, b2) &&
                            b2 != a1 &&
                            b2 != a2 &&
                            !g.content.edge(a1, b1) &&
                            !g.content.edge(a2, b1) &&
                            !g.content.edge(a1, b2) &&
                            !g.content.edge(a2, b2) {
                                return true
                            }
                    }
                }
            }
        };
        false
    })
}

pub fn main() {
    init_default_log();
    
    let n = 6;
    let basis = Basis::new(n);
    let edge: F = Graph::new(2, &[(0,1)]).into(); 
    // We optimize the average strong degree,
    // which is also the total number of edges (minus 1) in a strong clique
    let obj = ( strong_degree() * edge.project(n-2)).untype();
    let mut ineqs = vec![
        flags_are_nonnegative(basis),
        not_strong_clique(n).at_most(0)
    ];
    ineqs.push( one(n).at_most(1) );
    ineqs.append( &mut Degree::regularity(basis));

    let pb = Problem::<i64, _> {
        ineqs,
        cs: basis.all_cs(),
        obj: -obj,
    };

    let value = -pb.solve_csdp("strong_complete").expect("Failed running csdp");
    println!("Maximal number of edges: {:?} times Δ choose 2", value);
}
