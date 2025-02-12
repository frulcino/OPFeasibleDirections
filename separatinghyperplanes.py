#%%
import random
import itertools
from typing import List, Tuple, Optional
#%%
def generate_random_cuboid(n: int,
                           low: float = -10.0,
                           high: float = 10.0) -> List[Tuple[float, float]]:
    """
    Generate a random n-dimensional cuboid (box) by choosing
    n intervals [l_i, u_i], each with l_i < u_i.
    Returns a list of intervals [(l_1, u_1), ..., (l_n, u_n)].
    """
    intervals = []
    for _ in range(n):
        a = random.uniform(low, high)
        b = random.uniform(low, high)
        l_i, u_i = min(a, b), max(a, b)
        # Ensure we have a strict interval (l_i < u_i),
        # or handle the case if l_i == u_i as you prefer.
        if l_i == u_i:
            u_i = l_i + 1e-8  # tiny shift
        intervals.append((l_i, u_i))
    return intervals


def corners_of_cuboid(intervals: List[Tuple[float, float]]):
    """
    Return all 2^n vertices of the cuboid in n dimensions.
    Each interval [l_i, u_i] -> pick either l_i or u_i.
    """
    # For each dimension i, we have two choices: intervals[i][0] or intervals[i][1].
    # Use itertools.product to generate cartesian product.
    list_of_endpoints = [
        (intervals[i][0], intervals[i][1]) for i in range(len(intervals))
    ]
    for selection in itertools.product(*list_of_endpoints):
        yield selection  # a tuple (x_1, ..., x_n)


def product_of_coordinates(x: List[float]) -> float:
    """Compute p(x) = product of all coordinates of x."""
    prod_val = 1.0
    for val in x:
        prod_val *= val
    return prod_val


def define_pi_a(intervals: List[Tuple[float, float]],
                a: Tuple[float, ...]):
    """
    Given:
      - intervals: the n intervals [l_i, u_i]
      - a: a corner of the cuboid (a_1, ..., a_n)
    Return a callable function pi_a(x) implementing the unique affine map
    that agrees with p(x)=∏ x_i at 'a' and all corners adjacent to 'a'.
    
    In multilinear form, for x in R^n:
      pi_a(x) = p(a) + Σ [p(a^(i)) - p(a)] * (x_i - a_i) / (a^(i)_i - a_i)
    where a^(i) differs from a in the i-th coordinate (using the other endpoint).
    """
    n = len(intervals)
    # Precompute p(a):
    p_a = product_of_coordinates(a)
    
    # For each coordinate i, define a^(i) by flipping a_i to the other endpoint
    # in intervals[i], and compute p(a^(i)) - p(a).
    slopes_list = []
    for i in range(n):
        l_i, u_i = intervals[i]
        # the 'other' endpoint from a_i:
        if abs(a[i] - l_i) < 1e-15:
            # a_i == l_i, so the opposite is u_i
            a_i_opposite = u_i
        else:
            # a_i == u_i, so the opposite is l_i
            a_i_opposite = l_i
        
        # Construct a^(i)
        if a[i] != 0:
            slope_i = p_a/a[i]
        else:  #calculate the product of all but a_i
            a_i_list = list(a)
            a_i_list[i] = 1
            a_i_tuple = tuple(a_i_list)
            slope_i= product_of_coordinates(a_i_tuple)
        
        slopes_list.append(slope_i)
    
    def pi_of_x(x: List[float]) -> float:
        """
        Evaluate pi_a at the point x.
        """
        val = p_a
        # Add each linear term:
        for i in range(n):
            # (x[i] - a[i])
            contribution = slopes_list[i] * (x[i] - a[i])
            val += contribution
        return val
    
    return pi_of_x


def random_point_in_cuboid(intervals: List[Tuple[float, float]]) -> List[float]:
    """
    Generate a random point x in the cuboid, each x[i]
    uniform in [l_i, u_i].
    """
    return [
        random.uniform(l_i, u_i) for (l_i, u_i) in intervals
    ]


def is_separating_hyperplane(intervals: List[Tuple[float, float]],
                             a: Tuple[float, ...],
                             num_samples: int = 10000) -> bool:
    """
    Check whether pi_a is a separating hyperplane for p(x) = ∏ x_i
    on the cuboid 'intervals'.
    
    We do this *approximately* by sampling 'num_samples' random points.
    If we find any x in the cuboid with pi_a(x) > p(x), we return False.
    Otherwise, we return True.
    
    WARNING: This is a heuristic check. If no violation is found
    in random sampling, we cannot be 100% certain no violation exists.
    """
    pi_func = define_pi_a(intervals, a)

    upper = True
    lower = True
    for _ in range(num_samples):
        x = random_point_in_cuboid(intervals)
        p_x = product_of_coordinates(x)
        pi_x = pi_func(x)
        if pi_x > p_x:
            lower = False
        elif pi_x < p_x:
            upper = False
        if not upper and not lower:
            return False
    return True


def find_counterexample(intervals: List[Tuple[float, float]],
                       num_corners_check: int = 0,
                       num_samples: int = 10000) -> Optional[Tuple[float, ...]]:
    """
    Attempt to find a corner 'a' of the cuboid for which pi_a is NOT a
    separating hyperplane. We do so by random sampling in the interior.
    
    - If num_corners_check <= 0 or > 2^n, we check ALL corners.
    - If we find a corner 'a' whose pi_a is not separating,
      we return that corner immediately.
    - Otherwise return None if all tested corners appear separating in samples.
    """
    all_corners = list(corners_of_cuboid(intervals))
    n_corners = len(all_corners)
    
    if num_corners_check <= 0 or num_corners_check > n_corners:
        corners_to_check = all_corners
    else:
        random.shuffle(all_corners)
        corners_to_check = all_corners[:num_corners_check]
    
    for corner in corners_to_check:
        if not is_separating_hyperplane(intervals, corner, num_samples):
            return corner  # found a counterexample
    return None

#%%
# ------------------- DEMO USAGE -------------------
if __name__ == "__main__":

    #%% Example usage:
    # 1) Generate a random 3D cuboid
    cuboid_3d = generate_random_cuboid(n=5, low=-5, high=5)
    print("Random 3D cuboid intervals:", cuboid_3d)

    #%% 2) Pick a corner at random
    corners_3d = list(corners_of_cuboid(cuboid_3d))
    a_random = random.choice(corners_3d)
    print("Random corner a:", a_random)

    #%% 3) Build pi_a
    pi_a_func = define_pi_a(cuboid_3d, a_random)

    #%% 4) Check if pi_a is separating by random sampling
    is_sep = is_separating_hyperplane(cuboid_3d, a_random, num_samples=2000)
    print(f"Is pi_a a separating hyperplane? {is_sep}")

    #%% 5) Try to find a corner whose pi_a is NOT separating
    bad_corner = find_counterexample(cuboid_3d, num_corners_check=10, num_samples=2000)
    if bad_corner is not None:
        print("Found a corner whose pi_a is NOT separating:", bad_corner)
    else:
        print("No counterexample found among the tested corners.")
