import ddf.minim.*;
import ddf.minim.analysis.*;
import peasy.*;

VisualizationManager visManager;

void setup() {
    size(800, 800, P3D);
    visManager = new VisualizationManager(this);
    frameRate(60);
}

void draw() {
    visManager.update();
    visManager.render();
}

class VisualizationManager {
    private AudioProcessor audioProcessor;
    private Hypercube cube;
    private TimeHistoryManager timeHistory;
    private PeasyCam cam;
    private ParticleRenderer renderer;
    
    public VisualizationManager(PApplet parent) {
        cube = new Hypercube(256, 32);
        audioProcessor = new AudioProcessor(parent, cube);
        timeHistory = new TimeHistoryManager(cube.NumLayers, cube.GridSize);
        renderer = new ParticleRenderer(cube);
        
        initializeCamera(parent);
    }

    private void initializeCamera(PApplet parent) {
        cam = new PeasyCam(parent, 800);
        cam.setMinimumDistance(200);
        cam.setMaximumDistance(2000);
        parent.sphereDetail(6);
    }

    public void update() {
        audioProcessor.update();
        float currentFreq = audioProcessor.getCurrentFrequency();
        float stereoWidth = audioProcessor.getStereoWidth();
        timeHistory.update(currentFreq, stereoWidth);
    }

    public void render() {
        background(0);
        drawFrame();
        drawFPS();
        renderer.render(timeHistory, audioProcessor.getStereoWidth());
    }
    
    private void drawFrame() {
        stroke(255, 0, 255);
        noFill();
        box(cube.CubeSize);
    }
    
    private void drawFPS() {
        fill(255);
        textSize(16);
        textAlign(RIGHT, TOP);
        text("FPS: " + int(frameRate), cube.CubeSize / 2, cube.CubeSize / 2, cube.CubeSize / 2);
    }
}

class AudioProcessor {
    private Minim minim;
    private AudioOutput out;
    private FFT fftLeft, fftRight;
    private FrequencySweeper sweeper;
    
    final int SAMPLE_RATE = 44100;
    final int BUFFER_SIZE = 1024;
    
    public AudioProcessor(PApplet parent, Hypercube cube) {
        minim = new Minim(parent);
        out = minim.getLineOut(Minim.STEREO, BUFFER_SIZE);
        
        initializeFFT();
        sweeper = new FrequencySweeper(cube.YMin, cube.YMax);
        cube.InitSmoothSpectrum(fftLeft.specSize());
    }
    
    public void update() {
        sweeper.update(frameRate);
        float[] bufferLeft = new float[BUFFER_SIZE];
        float[] bufferRight = new float[BUFFER_SIZE];
        generateAudioBuffer(bufferLeft, bufferRight);
        processFFT(bufferLeft, bufferRight);
    }
    
    private void generateAudioBuffer(float[] bufferLeft, float[] bufferRight) {
        float phaseIncrement = TWO_PI * sweeper.getCurrentFrequency() / SAMPLE_RATE;
        float phase = sweeper.getPhase();
        float stereoWidth = sweeper.getStereoWidth();
        
        for (int i = 0; i < bufferLeft.length; i++) {
            float signal = 0.5 * sin(phase);
            
            bufferLeft[i] = signal * (1.0 - stereoWidth * 0.5);
            bufferRight[i] = signal * (1.0 + stereoWidth * 0.5);
            
            phase += phaseIncrement;
            if (phase > TWO_PI) phase -= TWO_PI;
        }
        
        sweeper.setPhase(phase);
    }
    
    private void initializeFFT() {
        fftLeft = new FFT(BUFFER_SIZE, SAMPLE_RATE);
        fftRight = new FFT(BUFFER_SIZE, SAMPLE_RATE);
        fftLeft.window(FFT.HAMMING);
        fftRight.window(FFT.HAMMING);
    }
    
    private void processFFT(float[] bufferLeft, float[] bufferRight) {
        fftLeft.forward(bufferLeft);
        fftRight.forward(bufferRight);
    }
    
    public float getCurrentFrequency() { return sweeper.getCurrentFrequency(); }
    public float getStereoWidth() { return sweeper.getStereoWidth(); }
}

class FrequencySweeper {
    private float currentFreq;
    private float minFreq;
    private float maxFreq;
    private float sweepDuration = 1.0;
    private float sweepTime = 0.0;
    private boolean sweepingUp = true;
    private float stereoWidth = 0.0;
    private float phase = 0.0;
    
    public FrequencySweeper(float minFreq, float maxFreq) {
        this.minFreq = minFreq;
        this.maxFreq = maxFreq;
        this.currentFreq = minFreq;
    }
    
    public void update(float frameRate) {
        float timeIncrement = 1.0 / frameRate;
        sweepTime += timeIncrement;
        
        if (sweepTime >= sweepDuration) {
            sweepTime = 0;
            sweepingUp = !sweepingUp;
        }
        
        float sweepProgress = sweepTime / sweepDuration;
        updateFrequencyAndStereo(sweepProgress);
    }
    
    private void updateFrequencyAndStereo(float progress) {
        if (sweepingUp) {
            currentFreq = minFreq * pow(maxFreq/minFreq, progress);
            stereoWidth = progress;
        } else {
            currentFreq = maxFreq * pow(minFreq/maxFreq, progress);
            stereoWidth = 1.0 - progress;
        }
    }
    
    public float getCurrentFrequency() { return currentFreq; }
    public float getStereoWidth() { return stereoWidth; }
    public float getPhase() { return phase; }
    public void setPhase(float phase) { this.phase = phase; }
}

class TimeHistoryManager {
    private TimeSlice[] timeHistory;
    private int currentSliceIndex = 0;
    private float lastUpdateTime = 0;
    private final int TIME_HISTORY = 1000; // milliseconds
    private ParticleGrid grid;
    private float currentStereoWidth = 0.0;
    
    public TimeHistoryManager(int numLayers, int gridSize) {
        timeHistory = new TimeSlice[numLayers];
        for (int i = 0; i < numLayers; i++) {
            timeHistory[i] = new TimeSlice(gridSize);
        }
        grid = new ParticleGrid(new Hypercube(gridSize * gridSize, numLayers));
    }
    
    public void update(float currentFreq, float stereoWidth) {
        currentStereoWidth = stereoWidth; // Store current stereo width
        float currentTime = millis();
        float deltaTime = currentTime - lastUpdateTime;
        
        if (deltaTime >= (TIME_HISTORY / timeHistory.length)) {
            currentSliceIndex = (currentSliceIndex + 1) % timeHistory.length;
            timeHistory[currentSliceIndex].clear();
            lastUpdateTime = currentTime;
        }
        
        // Calculate new particle sizes for current slice
        TimeSlice currentSlice = getCurrentSlice();
        updateParticleSizes(currentSlice, currentFreq, stereoWidth);
    }
    
    public float getCurrentStereoWidth() {
        return currentStereoWidth;
    }
    
    private void updateParticleSizes(TimeSlice slice, float currentFreq, float stereoWidth) {
        for (int row = 0; row < slice.getGridSize(); row++) {
            float rowFreq = grid.getRowFrequency(row);
            float freqDistance = abs(rowFreq - currentFreq);
            float relativeDistance = freqDistance / currentFreq;
            float freqFactor = (relativeDistance < 0.1) ? map(relativeDistance, 0, 0.1, 1.0, 0.0) : 0.0;
            
            for (int col = 0; col < slice.getGridSize(); col++) {
                float centerDistance = grid.getCenterDistance(col);
                
                // Calculate stereo factor
                float stereoFactor;
                if (stereoWidth < 0.1) {
                    stereoFactor = (1.0 - centerDistance) * map(stereoWidth, 0, 0.1, 1.0, 0.8);
                } else if (stereoWidth > 0.9) {
                    stereoFactor = centerDistance * map(stereoWidth, 0.9, 1.0, 0.8, 1.0);
                } else {
                    float monoComponent = (1.0 - centerDistance) * (1.0 - stereoWidth);
                    float stereoComponent = centerDistance * stereoWidth;
                    stereoFactor = max(monoComponent, stereoComponent);
                }
                
                // Calculate particle size
                float particleSize = ParticleGrid.BASE_SIZE;
                if (freqFactor > 0) {
                    particleSize += (ParticleGrid.MAX_SIZE - ParticleGrid.BASE_SIZE) * freqFactor * stereoFactor;
                }
                
                // Update the slice
                slice.setSizeAt(row, col, particleSize);
                slice.setStereoFactorAt(row, col, stereoFactor);
            }
        }
    }
    
    public TimeSlice getSlice(int layer) {
        int index = (timeHistory.length + currentSliceIndex - layer) % timeHistory.length;
        return timeHistory[index];
    }
    
    public TimeSlice getCurrentSlice() {
        return timeHistory[currentSliceIndex];
    }
    
    public int getCurrentIndex() {
        return currentSliceIndex;
    }
}

class ParticleRenderer {
    private Hypercube cube;
    private ParticleGrid grid;

    // Default value - higher = sharper peaks
    // Defines peak sharpness of stereo field
    private float fieldSeparation = 3.0;  
    
    public ParticleRenderer(Hypercube cube) {
        this.cube = cube;
        this.grid = new ParticleGrid(cube);
    }
    
    // Adjust how distinctly the peaks are separated
    public void setFieldSeparation(float separation) {
        this.fieldSeparation = separation;
    }
    
    public float getFieldSeparation() {
        return fieldSeparation;
    }
    
    private float calculateStereoFieldOpacity(float position, float stereoWidth) {
        // Convert position from 0-1 to -1 to 1 range (left to right)
        float stereoPosition = (position * 2.0) - 1.0;
        
        // Calculate where the signal should be strongest based on stereo width
        float peakPosition = stereoPosition < 0 ? -stereoWidth : stereoWidth;
        
        // Calculate distance from the peak
        float distance = abs(stereoPosition - peakPosition);
        
        // Use fieldSeparation to control the sharpness of the peaks
        return max(0, 1.0 - (distance * fieldSeparation));
    }
    
    public void render(TimeHistoryManager timeHistory, float stereoWidth) {
        translate(-(cube.CubeSize / 2), -(cube.CubeSize / 2), 0);
        noLights();
        noStroke();
        
        for (int layer = 0; layer < cube.NumLayers; layer++) {
            TimeSlice slice = timeHistory.getSlice(layer);
            float z = grid.getZPosition(layer);
            
            for (int row = 0; row < cube.GridSize; row++) {
                float y = grid.getYPosition(row);
                
                for (int col = 0; col < cube.GridSize; col++) {
                    float x = grid.getXPosition(col);
                    float size = slice.getSizeAt(row, col);
                    
                    if (size > ParticleGrid.BASE_SIZE) {
                        float xPos = (float)col / (cube.GridSize - 1);
                        float opacity = calculateStereoFieldOpacity(xPos, stereoWidth);
                        
                        if (opacity > 0.05) {
                            float alpha = 255 * opacity;
                            fill(0, 255, 255, alpha);
                            
                            pushMatrix();
                            translate(x, y, z);
                            sphere(size);
                            popMatrix();
                        }
                    }
                }
            }
        }
    }
}

class ParticleGrid {
    public static final float BASE_SIZE = 1.0;
    public static final float MAX_SIZE = 5.0;
    
    private float[] xPositions;
    private float[] yPositions;
    private float[] zPositions;
    private float[] rowFrequencies;
    private float[] stereoPositions;
    private float[] centerDistances;
    
    public ParticleGrid(Hypercube cube) {
        initializeGridPositions(cube);
    }
    
    private void initializeGridPositions(Hypercube cube) {
        xPositions = new float[cube.GridSize];
        yPositions = new float[cube.GridSize];
        zPositions = new float[cube.NumLayers];
        rowFrequencies = new float[cube.GridSize];
        stereoPositions = new float[cube.GridSize];
        centerDistances = new float[cube.GridSize];
        
        for (int i = 0; i < cube.GridSize; i++) {
            xPositions[i] = map(i, 0, cube.GridSize - 1, 0, cube.CubeSize);
            yPositions[i] = map(i, 0, cube.GridSize - 1, cube.CubeSize, 0);
            rowFrequencies[i] = cube.YMin * pow(cube.YMax/cube.YMin, (float)i/(cube.GridSize-1));
            stereoPositions[i] = (float)i / (cube.GridSize - 1);
            centerDistances[i] = abs(stereoPositions[i] - 0.5) * 2.0;
        }
        
        for (int i = 0; i < cube.NumLayers; i++) {
            zPositions[i] = map(i, 0, cube.NumLayers - 1, cube.CubeSize / 2, -cube.CubeSize / 2);
        }
    }
    
    public float getXPosition(int col) { return xPositions[col]; }
    public float getYPosition(int row) { return yPositions[row]; }
    public float getZPosition(int layer) { return zPositions[layer]; }
    public float getRowFrequency(int row) { return rowFrequencies[row]; }
    public float getCenterDistance(int col) { return centerDistances[col]; }
}

class TimeSlice {
    private float[][] sizes;
    private float[][] stereoFactors;
    
    public TimeSlice(int gridSize) {
        sizes = new float[gridSize][gridSize];
        stereoFactors = new float[gridSize][gridSize];
        clear();
    }
    
    public int getGridSize() {
        return sizes.length;
    }
    
    public void clear() {
        for (int i = 0; i < sizes.length; i++) {
            for (int j = 0; j < sizes[i].length; j++) {
                sizes[i][j] = 1.0;
                stereoFactors[i][j] = 0.0;
            }
        }
    }
    
    public void setSizeAt(int row, int col, float size) {
        sizes[row][col] = size;
    }
    
    public float getSizeAt(int row, int col) {
        return sizes[row][col];
    }
    
    public void setStereoFactorAt(int row, int col, float factor) {
        stereoFactors[row][col] = factor;
    }
    
    public float getStereoFactorAt(int row, int col) {
        return stereoFactors[row][col];
    }
}

class Hypercube {
    public static final float MIN_AMPLITUDE = 1e-6f;
    public static final float REF_AMPLITUDE = 1.0f;

    public int CubeSize = 400;
    public int NumParticles = 256;
    public int GridSize = (int) sqrt(NumParticles);
    public int NumLayers = 21;
    
    public float YMin = 50.0;    // Hz - Match sweep minimum
    public float YMax = 1000.0;  // Hz - Match sweep maximum
    public float[] smoothedSpectrum;
    
    public Hypercube(int numParticles, int numLayers) {
        NumParticles = numParticles;
        GridSize = (int) sqrt(numParticles);
        NumLayers = numLayers;
    }
    
    public void InitSmoothSpectrum(int specSize) {
        smoothedSpectrum = new float[specSize];
    }

    public float ToDb(float amplitude) {
        float absAmplitude = Math.abs(amplitude);
        float safeAmplitude = Math.max(absAmplitude, MIN_AMPLITUDE);
        return 20 * (float)Math.log10(safeAmplitude / REF_AMPLITUDE);
    }
}
