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
        s = "<h4>%s</h4>" % self.question
        s += "<form>"
        for ans in self.answers:
            s += f'<input type="radio" name="answer" value="{ans}">{ans}<br>'
        s += "</form>"

        answer = "La respuesta correcta es: <br>"
        answer += self.answers[self.correct_answer]
        if self.explanation is not None:
            answer += f"<br><i>" + self.explanation + "</i>"

        s += dedent(
            """
        <button title="Click to show/hide content" type="button"
        onclick="if(document.getElementById('{0}')
         .style.display=='none') {{document.getElementById('{0}')
         .style.display=''}}else{{document.getElementById('{0}')
         .style.display='none'}}">Show answer</button>
        <div id="{0}" style="display:none">{1}</div>"""
        )
        return s.format(secrets.token_urlsafe(20), answer)
