#!/usr/bin/env python3
# Classify rg scan lines (path:line:content) of the Alpöge archive into audit categories.
import sys, re, collections, os
SCAN, OUT = sys.argv[1], sys.argv[2]
KW = r'(syntax|macro|macro_rules|elab|elab_rules|notation|notation3|infix|infixl|infixr|prefix|postfix|declare_syntax_cat|binder_predicate|unif_hint)'
MOD = r'^\s*(@\[[^\]]*\]\s*)?((private|protected|scoped|local|noncomputable|partial|unsafe|meta)\s+)*'
r_syn = re.compile(MOD + KW + r'\b')
r_cmd = re.compile(MOD + r'(initialize|builtin_initialize|run_cmd|run_elab|run_meta|simproc|dsimproc|simproc_decl|register_simp_attr|register_option|register_builtin_option|irreducible_def|export|open private|add_decl_doc|declare_config_elab)\b')
r_attr = re.compile(r'@\[[^\]]*\b(macro|term_elab|command_elab|tactic|builtin_[a-z_]+|init|export|extern|implemented_by|csimp|command_parser|term_parser|tactic_parser|env_linter|delab|app_unexpander|incremental|computed_field|match_pattern|never_extract)\b')
r_api = re.compile(r'\b(addDecl|addDeclCore|addAndCompile|addDeclWithoutChecking|setEnv|modifyEnv|ofKernelEnv|unsafeCast|unsafeBaseIO|unsafeIO|unsafeEIO|ptrAddrUnsafe|importModules|withImportModules|setCheckedSync|addConstAsync|commitConst|commitCheckEnv|mkAuxTheorem|mkAuxLemma|mkAuxDefinition|replayConsts|evalConst|evalDecl|elabCommand|liftCommandElabM|liftTermElabM|ConstantVal|TheoremVal|DefinitionVal|thmDecl|defnDecl|axiomDecl|opaqueDecl|inductDecl|mkConst)\b|\b(Kernel|Environment|ConstantInfo|Declaration|IO|System)\.[A-Za-z]+')
r_bad = re.compile(r'\b(sorry|admit|sorryAx|axiom|native_decide|ofReduceBool|ofReduceNat|reduceBool|reduceNat|trustCompiler|opaque|implemented_by|extern|unsafe|skipKernelTC|Quot|Quotient)\b|_nested|debug\.')
r_opt = re.compile(r'set_option\s+([\w.]+)')
r_attrcmd = re.compile(r'^\s*attribute\s+\[([^\]]*)\]')
r_ind = re.compile(r'^\s*(mutual|inductive|class inductive)\b')
r_dec = re.compile(r'\b(decide|Decidable\.decide|decide_eq_true|of_decide_eq_true)\b')
r_wf = re.compile(r'\b(termination_by|decreasing_by|WellFounded\.fix|Acc\.rec|partial)\b')
r_num = re.compile(r'(?<![\w.])[0-9]{20,}|0x[0-9a-fA-F]{17,}')
r_k = re.compile(r'kexactB|kcopy|cpt_copy_thms|kc%|copyDec%')
cat = collections.defaultdict(list); cnt = collections.Counter(); opts = collections.Counter(); attrs = collections.Counter()
apis = collections.Counter(); bads = collections.Counter(); decf = collections.Counter(); wff = collections.Counter(); kf = collections.Counter()
def strip_comment(s):
    i = s.find('--')
    return s if i < 0 else s[:i]
with open(SCAN, encoding='utf-8', errors='replace') as f:
    for line in f:
        try:
            path, ln, txt = line.rstrip('\n').split(':', 2)
        except ValueError:
            continue
        loc = path[4:] + ':' + ln
        code = strip_comment(txt)
        if r_syn.match(code): cat['syntax'].append((loc, txt.strip()[:160]))
        if r_cmd.match(code): cat['cmd'].append((loc, txt.strip()[:160]))
        if r_attr.search(code): cat['attr'].append((loc, txt.strip()[:160]))
        for m in r_api.finditer(code):
            apis[m.group(0)] += 1
            if len(cat['api']) < 400000: cat['api'].append((loc, m.group(0), txt.strip()[:140]))
        for m in r_bad.finditer(code):
            bads[m.group(0)] += 1; cat['bad'].append((loc, m.group(0), txt.strip()[:140]))
        m = r_opt.search(code)
        if m: opts[m.group(1)] += 1; 
        if m and not m.group(1).startswith(('backward.isDefEq.respectTransparency','maxHeartbeats','synthInstance.maxHeartbeats','linter.','maxRecDepth')):
            cat['opt'].append((loc, m.group(1), txt.strip()[:140]))
        m = r_attrcmd.match(code)
        if m: attrs[m.group(1).split()[0] if m.group(1).split() else ''] += 1
        if r_ind.match(code): cat['ind'].append((loc, txt.strip()[:160]))
        if r_dec.search(code): cnt['decide_lines'] += 1; decf[path] += 1; cat['dec'].append((loc, txt.strip()[:140]))
        for m in r_wf.finditer(code): cnt['wf:' + m.group(1)] += 1; wff[path] += 1; cat['wf'].append((loc, m.group(1), txt.strip()[:140]))
        for m in r_num.finditer(code): cnt['bignum'] += 1; cat['num'].append((loc, m.group(0)[:60], txt.strip()[:140]))
        for m in r_k.finditer(code): kf[m.group(0)] += 1; cat['k'].append((loc, m.group(0), txt.strip()[:140]))
for k, v in cat.items():
    with open(os.path.join(OUT, 'cls_' + k + '.tsv'), 'w') as g:
        for row in v: g.write('\t'.join(row) + '\n')
print('syntax decls', len(cat['syntax']), '| cmds', len(cat['cmd']), '| attr', len(cat['attr']), '| inductive/mutual lines', len(cat['ind']))
print('cmd kinds', collections.Counter(re.sub(MOD, '', c[1]).split()[0] if re.sub(MOD,'',c[1]).split() else '' for c in cat['cmd']).most_common())
print('bad tokens', bads.most_common(30))
print('api top', apis.most_common(45))
print('set_option', opts.most_common(40))
print('attribute cmds', attrs.most_common(15))
print('decide lines', cnt['decide_lines'], 'files', len(decf), '| wf', {k: v for k, v in cnt.items() if k.startswith('wf:')}, 'files', len(wff), '| bignum', cnt['bignum'], '| k-tactics', kf.most_common())
