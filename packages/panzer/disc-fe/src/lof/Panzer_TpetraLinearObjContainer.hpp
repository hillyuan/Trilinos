// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_TpetraLinearObjContainer_hpp__
#define __Panzer_TpetraLinearObjContainer_hpp__

#include "PanzerDiscFE_config.hpp"

#include <map>

// Tpetra includes
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_applyDirichletBoundaryCondition.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"

#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_NodeType.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=panzer::TpetraNodeType>
class TpetraLinearObjContainer : public LinearObjContainer
                               , public ThyraObjContainer<ScalarT> {
   TpetraLinearObjContainer();

public:
   typedef LinearObjContainer::Members Members;

   typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
   typedef Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> CrsMatrixType;
   typedef Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> CrsGraphType;
   typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
   typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;
   typedef Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> ExportType;

   TpetraLinearObjContainer(const Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > & domain,
                            const Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > & range)
   {
      domainSpace = Thyra::createVectorSpace<ScalarT>(domain);
      rangeSpace = Thyra::createVectorSpace<ScalarT>(range);
   }

   virtual void initialize() 
   {
      if(get_x()!=Teuchos::null) get_x()->putScalar(0.0);
      if(get_dxdt()!=Teuchos::null) get_dxdt()->putScalar(0.0);
      if(get_f()!=Teuchos::null) get_f()->putScalar(0.0);
      if(get_A()!=Teuchos::null) {
        Teuchos::RCP<CrsMatrixType> mat = get_A(); 
        mat->resumeFill();
        mat->setAllToScalar(0.0);
        mat->fillComplete();
      }
   }

   //! Wipe out stored data.
   void clear()
   {
      set_x(Teuchos::null);
      set_dxdt(Teuchos::null);
      set_f(Teuchos::null);
      set_A(Teuchos::null);
   }

   inline void set_x(const Teuchos::RCP<VectorType> & in) { x = in; } 
   inline const Teuchos::RCP<VectorType> get_x() const { return x; }

   inline void set_dxdt(const Teuchos::RCP<VectorType> & in) { dxdt = in; } 
   inline const Teuchos::RCP<VectorType> get_dxdt() const { return dxdt; }

   inline void set_f(const Teuchos::RCP<VectorType> & in) { f = in; } 
   inline const Teuchos::RCP<VectorType> get_f() const { return f; }

   inline void set_A(const Teuchos::RCP<CrsMatrixType> & in) { A = in; } 
   inline const Teuchos::RCP<CrsMatrixType> get_A() const { return A; }

   void initializeMatrix(ScalarT value)
   {  
     A->resumeFill();
     A->setAllToScalar(value); 
     A->fillComplete();
   }

   virtual void set_x_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in) 
   { 
     if(in==Teuchos::null) { x = Teuchos::null; return; }

     Teuchos::RCP<const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > x_const 
         = TOE::getConstTpetraVector(in);
     x = Teuchos::rcp_const_cast<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(x_const); 
   } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_x_th() const 
   { return (x==Teuchos::null) ? Teuchos::null : Thyra::createVector(x,domainSpace); }

   virtual void set_dxdt_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { 
     if(in==Teuchos::null) { dxdt = Teuchos::null; return; }

     Teuchos::RCP<const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > dxdt_const 
         = TOE::getConstTpetraVector(in);
     dxdt = Teuchos::rcp_const_cast<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(dxdt_const); 
   } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_dxdt_th() const 
   { return (dxdt==Teuchos::null) ? Teuchos::null : Thyra::createVector(dxdt,domainSpace); }

   virtual void set_f_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { f = (in==Teuchos::null) ? Teuchos::null : TOE::getTpetraVector(in); } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_f_th() const 
   { return (f==Teuchos::null) ? Teuchos::null : Thyra::createVector(f,rangeSpace); }

   virtual void set_A_th(const Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > & in) 
   { A = (in==Teuchos::null) ? Teuchos::null : Teuchos::rcp_dynamic_cast<CrsMatrixType>(TOE::getTpetraOperator(in),true); }
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > get_A_th() const
   { return (A==Teuchos::null) ? Teuchos::null : Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(A,rangeSpace,domainSpace); }

   void applyDirichletBoundaryCondition( const std::map< panzer::LocalOrdinal, double >& indx ) override
   {
	   /*using device_type = typename CrsMatrixType::device_type;
      using execution_space = typename CrsMatrixType::execution_space;
      using range_type = Kokkos::RangePolicy<execution_space, LocalOrdinalT>;

      std::vector<panzer::LocalOrdinal> lids;
      for( auto itr: indx ) lids.emplace_back(itr.first);
	  
	   const LocalOrdinalT lclNumRows = indx.size();
	   Kokkos::View<typename CrsMatrixType::local_ordinal_type*, device_type> lclRowInds ("lclRowInds", lclNumRows);
	   Kokkos::parallel_for
      ("Fill lclRowInds",
         range_type (0, lclNumRows),
         KOKKOS_LAMBDA (const LocalOrdinalT lclRow) {
	        lclRowInds(lclRow) = lids[lclRow];
         }
      );
	   
	   Tpetra::applyDirichletBoundaryConditionToLocalMatrixRows(*A, lclRowInds);*/

    //  ContainerType & t_ghosted = Teuchos::dyn_cast<ContainerType>(GhostedContainer);
    //  Teuchos::RCP<CrsMatrixType> A = t_ghosted.get_A();
    //  VectorType f = *( t_ghosted.get_f() );
    //  Teuchos::ArrayRCP<ScalarT> f_1dview = f.get1dViewNonConst();

      for( auto itr: indx )
      {
      //   std::cout << itr.first << "," << itr.second << std::endl;
	      std::size_t numEntries = 0;
         std::size_t sz = A->getNumEntriesInLocalRow(itr.first);
         Teuchos::Array<LocalOrdinalT> indices(sz);
         Teuchos::Array<ScalarT> Entries(sz);
         A->getLocalRowCopy(itr.first,indices,Entries,numEntries);
	      for (std::size_t i=0; i<sz; i++) {
		      if( indices[i]==itr.first )
			      Entries[i] = 1.0;
		      else
		  	      Entries[i] = 0.0;
	      }  
         A->replaceLocalValues(itr.first,indices,Entries);
      //   f_1dview[itr.first] = 0.0;
      }
   }
   
   void evalDirichletResidual( const std::map< panzer::LocalOrdinal, double >& indx ) override
   {
   /*   using device_type = typename CrsMatrixType::device_type;
      using execution_space = typename CrsMatrixType::execution_space;
      using range_type = Kokkos::RangePolicy<execution_space, LocalOrdinalT>;
      
      const LocalOrdinalT lclNumRows = indx.size();
	   Kokkos::View<typename CrsMatrixType::local_ordinal_type*, device_type> lclRowInds ("lclRowInds", lclNumRows);
	   Kokkos::parallel_for
      ("Fill lclRowInds",
         range_type (0, lclNumRows),
         KOKKOS_LAMBDA (const LocalOrdinalT lclRow) {
	        f->replaceLocalValue(lclRow, 0.0);
         }
      );*/

	  Teuchos::ArrayRCP<const ScalarT> x_1dview = x->get1dView();
      for( auto itr: indx )
      {
	  	double a = x_1dview[itr.first] - itr.second;
         f->replaceLocalValue(itr.first, a);
      }
   }
   
   void writeMatrixMarket(const std::string& filename) const override
   {
	   Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>>
	   		::writeSparseFile(filename, *A);
   }
    
private:
   typedef Thyra::TpetraOperatorVectorExtraction<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> TOE;

   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > domainSpace;
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > rangeSpace;

   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > x, dxdt, f;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > A;
};

}

#endif
