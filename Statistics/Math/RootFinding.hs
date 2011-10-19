{-# LANGUAGE BangPatterns, DeriveDataTypeable #-}

-- |
-- Module    : Statistics.Math.RootFinding
-- Copyright : (c) 2011 Bryan O'Sullivan
-- License   : BSD3
--
-- Maintainer  : bos@serpentine.com
-- Stability   : experimental
-- Portability : portable
--
-- Haskell functions for finding the roots of mathematical functions.

module Statistics.Math.RootFinding
    (
      Root(..)
    , fromRoot
    , ridders
    -- * References
    -- $references
    ) where

import Control.Monad (MonadPlus(..), ap)
import Control.Applicative
import Data.Typeable (Typeable)

-- | The result of searching for a root of a mathematical function.
data Root a = NotBracketed
            -- ^ The function does not have opposite signs when
            -- evaluated at the lower and upper bounds of the search.
            | SearchFailed
            -- ^ The search failed to converge to within the given
            -- error tolerance after the given number of iterations.
            | Root a
            -- ^ A root was successfully found.
              deriving (Eq, Read, Show, Typeable)

instance Functor Root where
    fmap _ NotBracketed = NotBracketed
    fmap _ SearchFailed = SearchFailed
    fmap f (Root a)     = Root (f a)

instance Monad Root where
    NotBracketed >>= _ = NotBracketed
    SearchFailed >>= _ = SearchFailed
    Root a       >>= m = m a

    return = Root

instance MonadPlus Root where
    mzero = SearchFailed

    r@(Root _) `mplus` _ = r
    _          `mplus` p = p

instance Applicative Root where
    pure  = Root
    (<*>) = ap

instance Alternative Root where
    empty = SearchFailed

    r@(Root _) <|> _ = r
    _          <|> p = p

-- | Returns either the result of a search for a root, or the default
-- value if the search failed.
fromRoot :: a                   -- ^ Default value.
         -> Root a              -- ^ Result of search for a root.
         -> a
fromRoot _ (Root a) = a
fromRoot a _        = a

-- | Use the method of Ridders to compute a root of a function.
--
-- The function must have opposite signs when evaluated at the lower
-- and upper bounds of the search (i.e. the root must be bracketed).
ridders :: Double               -- ^ Error tolerance.
        -> (Double,Double)      -- ^ Lower and upper bounds for the search.
        -> (Double -> Double)   -- ^ Function to find the roots of.
        -> Root Double
ridders tol (lo,hi) f
    | flo ~= 0    = Root lo
    | fhi ~= 0    = Root hi
    | flo*fhi > 0 = NotBracketed -- root is not bracketed
    | otherwise   = go lo flo hi fhi 0
  where
    go !a !fa !b !fb !i
        | fn ~= 0 || abs (b-a) < tol = Root n
        | i >= (50 :: Int)           = SearchFailed
        | fn*fm < 0 = go n fn m fm (i+1)
        | fn*fa < 0 = go a fa n fn (i+1)
        | otherwise = go n fn b fb (i+1)
      where
        dm  = (b - a) * 0.5
        !m   = a + dm
        !fm  = f m
        !dn  = signum (fb - fa) * dm * fm / sqrt(fm*fm - fa*fb)
        !n   = m - signum dn * min (abs dn) (abs dm - 0.5 * tol)
        !fn  = f n
    !flo = f lo
    !fhi = f hi
    a ~= b = abs (a-b) <= abs (a*tol)

-- $references
--
-- * Ridders, C.F.J. (1979) A new algorithm for computing a single
--   root of a real continuous function.
--   /IEEE Transactions on Circuits and Systems/ 26:979&#8211;980.