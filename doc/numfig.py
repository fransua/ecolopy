from docutils.nodes import figure, caption, Text, reference, raw, SkipNode, Element
from sphinx.roles import XRefRole

def setup(app):
    app.add_config_value('number_figures', True, True)
    app.add_config_value('figure_caption_prefix', "Figure", True)

    app.add_node(page_ref,
                 text=(skip_page_ref, depart_page_ref),
                 html=(skip_page_ref, depart_page_ref),
                 latex=(visit_page_ref, depart_page_ref))

    app.add_role('page', XRefRole(nodeclass=page_ref))

    app.add_node(num_ref,
                 text=(visit_num_ref, depart_num_ref),
                 html=(visit_num_ref, depart_num_ref),
                 latex=(latex_visit_num_ref, latex_depart_num_ref))

    app.add_role('num', XRefRole(nodeclass=num_ref))

    app.connect('doctree-resolved', number_figure_nodes)

# Element classes
class page_ref(reference):
    pass

class num_ref(reference):
    pass

def skip_page_ref(self, node):
    raise SkipNode

def visit_page_ref(self, node):
    self.body.append("\\pageref{%s:%s}" % (node['refdoc'], node['reftarget']))
    raise SkipNode

def depart_page_ref(self, node):
    pass

def visit_num_ref(self, node):
    pass

def depart_num_ref(self, node):
    pass

def latex_visit_num_ref(self, node):
    fields = node['reftarget'].split('#')
    if len(fields) > 1:
        label, target = fields
        ref_link = '%s:%s' % (node['refdoc'], target)
        latex = "\\hyperref[%s]{%s \\ref*{%s}}" % (ref_link, label, ref_link)
        self.body.append(latex)
    else:
        self.body.append('\\ref{%s:%s}' % (node['refdoc'], fields[0]))
        
    raise SkipNode

def latex_depart_num_ref(self, node):
    pass

def number_figure_nodes(app, doctree, docname):
    env = app.builder.env
    # first generate figure numbers for each figure
    i = 1
    figids = {}
    for figure_info in doctree.traverse(figure):
        if app.builder.name != 'latex' and app.config.number_figures:
            for cap in figure_info.traverse(caption):
                cap[0] = Text("%s %d: %s" % (app.config.figure_caption_prefix, i, cap[0]))

        for id in figure_info['ids']:
            figids[id] = i
        i += 1

    # replace numfig nodes with links
    if app.builder.name != 'latex':
        for ref_info in doctree.traverse(num_ref):
            if '#' in ref_info['reftarget']:
                label, target = ref_info['reftarget'].split('#')
                labelfmt = label + " %d"
            else:
                labelfmt = '%d'
                target = ref_info['reftarget']

            if app.builder.name == 'html':
                link = ref_info['refdoc']+'.html#'+target
                html = '<a href="%s">%s</a>' % (link, labelfmt %(figids[target]))
                ref_info.replace_self(raw(html, html, format='html'))
            else:
                ref_info.replace_self(Text(labelfmt % (figids[target])))
