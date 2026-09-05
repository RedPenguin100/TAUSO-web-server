import os
import io
import zipfile
import base64
import logging
from html import escape
from dotenv import load_dotenv
import brevo_python
from brevo_python.rest import ApiException

load_dotenv(".env.local")

BREVO_API_KEY = os.getenv("BREVO_API_KEY")
FROM_EMAIL = os.getenv("FROM_EMAIL")
FROM_NAME = os.getenv("FROM_NAME")

# Keep a runaway exception message from bloating the failure email.
REASON_MAX_CHARS = 300

logger = logging.getLogger(__name__)

def send_mail(to_email: str, subject: str, html_content: str, attachments: list = None):
    """
    Synchronous function to send emails via Brevo API with optional zip attachments.
    attachments: list of tuples (filename, bytes_content)
    """
    if not BREVO_API_KEY:
        logger.error("BREVO_API_KEY is missing!")
        return

    configuration = brevo_python.Configuration()
    configuration.api_key['api-key'] = BREVO_API_KEY
    api_instance = brevo_python.TransactionalEmailsApi(brevo_python.ApiClient(configuration))

    email_attachments = []

    # Zip attachments if any are provided (matching your FastAPI logic)
    if attachments:
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for filename, content in attachments:
                zip_file.writestr(filename, content)

        zip_bytes = zip_buffer.getvalue()
        base64_content = base64.b64encode(zip_bytes).decode('utf-8')

        email_attachments.append(
            brevo_python.SendSmtpEmailAttachment(
                content=base64_content,
                name="tauso_analysis_results.zip"
            )
        )

    send_smtp_email = brevo_python.SendSmtpEmail(
        sender={"email": FROM_EMAIL, "name": FROM_NAME},
        to=[{"email": to_email}],
        subject=subject,
        html_content=html_content,
        attachment=email_attachments if email_attachments else None
    )

    try:
        api_instance.send_transac_email(send_smtp_email)
        logger.info(f"Email successfully sent to {to_email}")
    except ApiException as e:
        logger.error(f"Exception when calling Brevo API: {e}")

def send_processing_started(to_email: str, source_info: str):
    html = f"""
    <div style="font-family: Arial, sans-serif;">
      <h2>Pipeline Initialized</h2>
      <p>We've started processing your sequence via TAUSO.</p>
      <p><strong>Source:</strong> {escape(source_info)}</p>
      <p>You will receive another email once the ASO efficacy prediction is complete.</p>
    </div>
    """
    send_mail(to_email, "TAUSO Pipeline Started", html)

def send_processing_completed(to_email: str, source_info: str, results_files: list):
    html = f"""
    <div style="font-family: Arial, sans-serif;">
      <h2>Analysis Complete! 🎉</h2>
      <p>Your TAUSO ASO efficacy predictions are ready.</p>
      <p><strong>Source:</strong> {escape(source_info)}</p>
      <p>Please find your results attached in the zip file.</p>
    </div>
    """
    send_mail(to_email, "TAUSO Analysis Complete", html, attachments=results_files)

def send_processing_failed(to_email: str, source_info: str, reason: str):
    """Tell the submitter their run produced no results. `reason` is summarised rather than a
    traceback, which would expose server internals to whoever submitted the job."""
    html = f"""
    <div style="font-family: Arial, sans-serif;">
      <h2>Analysis Failed</h2>
      <p>Your TAUSO run stopped before producing results.</p>
      <p><strong>Source:</strong> {escape(source_info)}</p>
      <p><strong>Reason:</strong> {escape(reason[:REASON_MAX_CHARS])}</p>
      <p>The failure has been logged. Re-submitting is worth trying if the cause was transient.</p>
    </div>
    """
    send_mail(to_email, "TAUSO Analysis Failed", html)