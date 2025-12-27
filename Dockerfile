FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy requirements file
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY calculate_volumes.py .
COPY api_service.py .
COPY index.html .

# Expose the API port
EXPOSE 5000

# Set environment variable for Flask
ENV FLASK_APP=api_service.py
ENV FLASK_ENV=production

# Run the API service
CMD ["python", "api_service.py"]

