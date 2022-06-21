from textwrap import dedent
import secrets
class MultipleChoice:
    def __init__(
        self,
        question,
        answers,
        correct_answer,
        explanation=None,
    ):
        """
        Multiple choice question html component
        Parameters:
        -----------
        question : string
        answers : list of strings
        correct_answers : int
        """
        self.question = question
        self.answers = answers
        self.correct_answer = correct_answer
        self.explanation = explanation

    def _repr_html_(self):
        s = '<div class="note admonition">'
        s+= '<p class="admonition-title">'
        #s+= '::before'
        s+= 'Pregunta'
        s+= '</p>'
        s += '<p><span class="math notranslate nohighlight">%s</span></p>' % self.question
        s += '<form>'
        for ans in self.answers:
            s += f'<input type="radio" name="answer" value="{ans}"><span class="math notranslate nohighlight">{ans}</span><br>'
        s += '</form>'
        answer = 'La respuesta correcta es: <br>'
        answer += f"<b>"+self.answers[self.correct_answer]+f"</b>"
        if self.explanation is not None:
            answer += f'<br><span class="math notranslate nohighlight">' + self.explanation + f'</span>'
        s += dedent(
"""
<details class="toggle-details">
<summary class="toggle-details__summary">
<svg xmlns="http://www.w3.org/2000/svg" class="tb-icon toggle-chevron-right" width="44" height="44" viewBox="0 0 24 24" stroke-width="1.5" stroke="#000000" fill="none" stroke-linecap="round" stroke-linejoin="round">
<path stroke="none" d="M0 0h24v24H0z" fill="none"></path>
<polyline points="9 6 15 12 9 18"></polyline>
</svg>
<span class="toggle-details__summary-text">Mostrar respuesta</span>
</summary>
<div class="cell_output docutils container">
<div><span class="math notranslate nohighlight">{0}</span></div>
</div></details>
"""
        )
        #s += dedent(
        #    """
        #<button title="Click to show/hide content" type="button"
        #onclick="if(document.getElementById('{0}')
        # .style.display=='none') {{document.getElementById('{0}')
        # .style.display=''}}else{{document.getElementById('{0}')
        # .style.display='none'}}">Mostrar la respuesta</button>
        #<div id="{0}" style="display:none">{1}</div>"""
        #)
        #s+='::after'
        s+='</div)'
        #return s.format(secrets.token_urlsafe(20), answer)
        return s.format(answer)
