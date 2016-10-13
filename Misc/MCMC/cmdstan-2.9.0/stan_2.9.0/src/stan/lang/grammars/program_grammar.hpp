#ifndef STAN_LANG_GRAMMARS_PROGRAM_GRAMMAR_HPP
#define STAN_LANG_GRAMMARS_PROGRAM_GRAMMAR_HPP

#include <boost/spirit/include/qi.hpp>

#include <stan/lang/ast.hpp>
#include <stan/lang/grammars/whitespace_grammar.hpp>
#include <stan/lang/grammars/expression_grammar.hpp>
#include <stan/lang/grammars/var_decls_grammar.hpp>
#include <stan/lang/grammars/statement_grammar.hpp>
#include <stan/lang/grammars/functions_grammar.hpp>

#include <string>
#include <sstream>
#include <utility>
#include <vector>

namespace stan {
  namespace lang {

    template <typename Iterator>
    struct program_grammar
      : boost::spirit::qi::grammar<Iterator,
                                   program(),
                                   whitespace_grammar<Iterator> > {
      program_grammar(const std::string& model_name);

      // global info for parses
      std::string model_name_;
      variable_map var_map_;
      std::stringstream error_msgs_;

      // grammars
      expression_grammar<Iterator> expression_g;
      var_decls_grammar<Iterator> var_decls_g;
      statement_grammar<Iterator> statement_g;
      functions_grammar<Iterator> functions_g;

      // rules

      boost::spirit::qi::rule<Iterator,
                              std::vector<var_decl>(),
                              whitespace_grammar<Iterator> >
      data_var_decls_r;

      boost::spirit::qi::rule<Iterator,
                              std::pair<std::vector<var_decl>,
                                        std::vector<statement> >(),
                              whitespace_grammar<Iterator> >
      derived_data_var_decls_r;

      boost::spirit::qi::rule<Iterator,
                              std::pair<std::vector<var_decl>,
                                        std::vector<statement> >(),
                              whitespace_grammar<Iterator> >
      derived_var_decls_r;

      boost::spirit::qi::rule<Iterator,
                              std::pair<std::vector<var_decl>,
                                        std::vector<statement> >(),
                              whitespace_grammar<Iterator> >
      generated_var_decls_r;

      boost::spirit::qi::rule<Iterator,
                              statement(),
                              whitespace_grammar<Iterator> >
      model_r;

      boost::spirit::qi::rule<Iterator,
                              std::vector<var_decl>(),
                              whitespace_grammar<Iterator> >
      param_var_decls_r;


      boost::spirit::qi::rule<Iterator,
                              program(),
                              whitespace_grammar<Iterator> >
      program_r;

      boost::spirit::qi::rule<Iterator,
                              boost::spirit::qi::unused_type,
                              whitespace_grammar<Iterator> >
      end_var_decls_r;

      boost::spirit::qi::rule<Iterator,
                              boost::spirit::qi::unused_type,
                              whitespace_grammar<Iterator> >
      end_var_decls_statements_r;

      boost::spirit::qi::rule<Iterator,
                              boost::spirit::qi::unused_type,
                              whitespace_grammar<Iterator> >
      end_var_definitions_r;
    };

  }
}
#endif
