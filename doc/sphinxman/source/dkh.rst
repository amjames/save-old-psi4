
.. include:: autodoc_abbr_options_c.rst

.. index:: DKH
.. _`sec:DKH`:

Interface to DKH by A. Wolf, M. Reiher, and B. A. Hess
======================================================

.. codeauthor:: Justin M. Turney
.. sectionauthor:: Justin M. Turney

*Module:*  :source:`DKH <src/lib/libmints>`

.. _`sec:dkhinput`:

Input
~~~~~

For all electron calculations one can use the Douglas-Kroll-Hess (DKH)
Hamiltonian to take into account scalar relativistic effects.

Minimal input for DKH single-point computation looks like this::

    molecule {
    Mg
    }

    set basis aug-cc-pvdz-dk
    set relativistic dkh

    energy('scf')

By default a 2nd-order DKH calculation is performed. To change the default
order use the |globals__dkh_order| option. The version of the code found in
|Psifour| is capable of up to 4th-order DKH calculations.

Keywords
~~~~~~~~

.. include:: autodir_options_c/globals__relativistic.rst
.. include:: autodir_options_c/globals__dkh_order.rst

Reference
~~~~~~~~~

When using this code please make reference to the appropriate following paper:

* "The Generalized Douglas-Kroll Transformation," A. Wolf,
  M. Reiher, and B. A. Hess, *J. Chem. Phys.* **117**, 9215 (2002).
  (doi: `10.1063/1.1515314 <http://dx.doi.org/10.1063/1.1515314>`_)

