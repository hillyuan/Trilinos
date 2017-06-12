#ifndef TEUCHOS_YAML_HPP
#define TEUCHOS_YAML_HPP

#include <Teuchos_Language.hpp>
#include <Teuchos_ReaderTables.hpp>

namespace Teuchos {
namespace yaml {

enum {
  PROD_DOC,
  PROD_PROLOG,
  PROD_NO_DIRECT,
  PROD_ONE_DIRECT,
  PROD_DIRECT,
  PROD_EPILOG,
  PROD_DOC_END,
  PROD_TOP_BMAP,
  PROD_TOP_BSEQ,
  PROD_TOP_BLOCK,
  PROD_BMAP,
  PROD_BSEQ,
  PROD_BMAP_ONE_ITEM,
  PROD_BMAP_ITEMS,
  PROD_BSEQ_ITEM,
  PROD_BSEQ_ITEMS,
  PROD_BSEQ_SCALAR,
  PROD_BMAP_ITEM,
  PROD_BMAP_SCALAR,
  PROD_BMAP_BLOCK,
  PROD_BMAP_FLOW,
  PROD_FSEQ_EMPTY,
  PROD_FMAP_EMPTY,
  PROD_FSEQ,
  PROD_FMAP,
  PROD_FSEQ_ITEM,
  PROD_FSEQ_ITEMS,
  PROD_FSEQ_SCALAR,
  PROD_FSEQ_FLOW,
  PROD_FMAP_ONE_ITEM,
  PROD_FMAP_ITEMS,
  PROD_FMAP_ITEM,
  PROD_FMAP_SCALAR,
  PROD_FMAP_FLOW,
  PROD_NO_SPACE,
  PROD_SPACE,
  PROD_NO_COMMENTS,
  PROD_COMMENTS,
  PROD_NO_EQDENT,
  PROD_EQDENT,
  PROD_RAW,
  PROD_DQUOTED,
  PROD_SQUOTED
};

enum { NPRODS = PROD_SQUOTED + 1 };

enum {
  TOK_NODENT,
  TOK_INDENT,
  TOK_DEDENT,
  TOK_EQDENT,
  TOK_SPACES,
  TOK_COMMENT,
  TOK_COLON,
  TOK_DOC_START,
  TOK_DOC_END,
  TOK_DIRECT,
  TOK_RAW,
  TOK_DQUOTED,
  TOK_SQUOTED,
  TOK_BSEQ,
  TOK_FSEP,
  TOK_LSQUARE,
  TOK_RSQUARE,
  TOK_LCURLY,
  TOK_RCURLY
};

enum { NTOKS = TOK_RCURLY + 1 };

Language make_language();
LanguagePtr ask_language();
ReaderTablesPtr ask_reader_tables();

}  // end namespace yaml
}  // end namespace Teuchos

#endif
