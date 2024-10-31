import ddf.minim.*;
import ddf.minim.analysis.*;
import peasy.*;
import processing.sound.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.ArrayList;

final int CUBE_SIZE = 400;
final int PARTICLES_PER_LAYER = 512;
final int LAYERS = 64;
final float Y_MIN = 50.0f;  // Hz
final float Y_MAX = 1000.0f; // Hz
final int SAMPLE_RATE = 44100;
// Use a larger buffer size for better frequency resolution
// 4096 gives us ~10.8Hz resolution (44100/4096) which is
// helpful for resolution even when the yMin is 50Hz
final int BUFFER_SIZE = 4096;

VisualizationManager visManager;

void setup() {
    // Can't use static vars here for some reason
    size(2000, 1600, P3D);
    visManager = new VisualizationManager(this);
    frameRate(60);
}

void draw() {
    visManager.update();
    visManager.render();
}

class VisualizationManager {
    //private TestSignalAnalyzer audioAnalyzer;
    private AudioFileAnalyzer audioAnalyzer;

    private PeasyCam cam;
    private Hypercube cube;
    private TimeHistoryManager timeHistory;
    private ParticleRenderer renderer;
    
    public VisualizationManager(PApplet parent) {
        cube = new Hypercube(CUBE_SIZE, PARTICLES_PER_LAYER, LAYERS, Y_MIN, Y_MAX);
        //audioAnalyzer = new TestSignalAnalyzer(parent, cube.yMin, cube.yMax, SAMPLE_RATE, BUFFER_SIZE);
        audioAnalyzer = new AudioFileAnalyzer(parent, cube.yMin, cube.yMax, SAMPLE_RATE, BUFFER_SIZE);
        audioAnalyzer.loadFile("test.wav");
        audioAnalyzer.play(true);

        cube.initSmoothSpectrum(audioAnalyzer.getSpecSize());

        timeHistory = new TimeHistoryManager(cube);
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
        audioAnalyzer.update();
        float currentFreq = audioAnalyzer.getCurrentFrequency();
        float stereoWidth = audioAnalyzer.getStereoWidth();
        timeHistory.update(currentFreq, stereoWidth);
    }

    public void render() {
        background(0);
        
        // Draw cube frame at origin
        stroke(33, 33, 33);
        strokeWeight(2);
        noFill();
        box(cube.cubeSize);
        
        // Render particles
        renderer.render(timeHistory, audioAnalyzer.getStereoWidth());
        
        drawFPS();
    }
    
    private void drawFrame() {
        stroke(255, 0, 255);
        noFill();
        box(cube.cubeSize);
    }
    
    private void drawFPS() {
        fill(255);
        textSize(16);
        textAlign(RIGHT, TOP);
        text("FPS: " + int(frameRate), cube.cubeSize / 2, cube.cubeSize / 2, cube.cubeSize / 2);
    }
}

class Hypercube {
    public static final float MIN_AMPLITUDE = 1e-6f;
    public static final float REF_AMPLITUDE = 1.0f;

    public static final float POINT_SIZE_MIN = 2.0;
    public static final float POINT_SIZE_MAX = 10.0;

    public int cubeSize = 400;
    public int numParticles = 256;
    public int gridSize = (int) sqrt(numParticles);
    public int numLayers = 21;
    
    public float yMin = 50.0;    // Hz - Match sweep minimum
    public float yMax = 1000.0;  // Hz - Match sweep maximum
    public float[] smoothedSpectrum;
    
    public Hypercube(int cubeSize, int numParticles, int numLayers, float yMin, float yMax) {
        this.numParticles = numParticles;
        this.gridSize = (int) sqrt(numParticles);
        this.numLayers = numLayers;
        this.yMin = yMin;
        this.yMax = yMax;
    }
    
    public void initSmoothSpectrum(int specSize) {
        smoothedSpectrum = new float[specSize];
    }

    public final float toDb(float amplitude) {
        float absAmplitude = Math.abs(amplitude);
        float safeAmplitude = Math.max(absAmplitude, MIN_AMPLITUDE);
        return 20 * (float)Math.log10(safeAmplitude / REF_AMPLITUDE);
    }
}

abstract class BaseAudioAnalyzer {
    protected float yMin;
    protected float yMax;
    protected int sampleRate;
    protected int bufferSize;

    protected float[] bufferLeft;
    protected float[] bufferRight;

    protected Minim minim;
    protected ddf.minim.analysis.FFT fftLeft, fftRight;

    protected float currentFrequency = 0.0f;
    protected float currentStereoWidth = 0.0f;
    private float prevStereoWidth = 0.0f;
    private float prevFrequency = 0.0f;
    private static final float FREQ_SMOOTH = 0.3f;
    
    public BaseAudioAnalyzer(PApplet parent, float yMin, float yMax, int sampleRate, int bufferSize) {
        this.yMin = yMin;
        this.yMax = yMax;
        this.sampleRate = sampleRate;

        this.bufferSize = bufferSize;
        
        bufferLeft = new float[this.bufferSize];
        bufferRight = new float[this.bufferSize];

        minim = new Minim(parent);
        initializeFFT();
        
        float binResolution = (float)sampleRate / this.bufferSize;
        System.out.println("=== FFT Parameters ===");
        System.out.println("Y-axis range: " + yMin + " to " + yMax + " Hz");
        System.out.println("Sample rate: " + sampleRate + " Hz");
        System.out.println("Buffer size: " + this.bufferSize);
        System.out.println("FFT spec size: " + fftLeft.specSize());
        System.out.println("Frequency resolution: " + binResolution + " Hz/bin");
        System.out.println("Minimum detectable frequency: " + binResolution + " Hz");
        System.out.println("Time window: " + (1000f * this.bufferSize / sampleRate) + " ms");
        System.out.println("====================");
    }
    
    public abstract void play(boolean loop);
    public abstract void pause();
    public abstract void update();
    
    public int getSpecSize() { return fftLeft.specSize(); }
    public float getCurrentFrequency() { return currentFrequency; }
    public float getStereoWidth() { return currentStereoWidth; }
    
    protected void initializeFFT() {
        fftLeft = new ddf.minim.analysis.FFT(bufferSize, sampleRate);
        fftRight = new ddf.minim.analysis.FFT(bufferSize, sampleRate);
        
        // Use Blackman-Harris window for better frequency separation
        fftLeft.window(ddf.minim.analysis.FFT.BLACKMAN);
        fftRight.window(ddf.minim.analysis.FFT.BLACKMAN);
        
        // Increase average bands for smoother analysis
        fftLeft.logAverages(22, 3);  // Start at 22Hz with 3 bands per octave
        fftRight.logAverages(22, 3);
    }
    
    protected void processFFT(float[] bufferLeft, float[] bufferRight) {
        fftLeft.forward(bufferLeft);
        fftRight.forward(bufferRight);
        
        findBasicFrequency();
        currentStereoWidth = calculateStereoWidth();
    }
    
    private void findBasicFrequency() {
        float maxAmp = 0;
        int maxBin = 0;
        
        // Start from bin 1 to avoid DC component
        for (int i = 1; i < fftLeft.specSize(); i++) {
            float leftAmp = fftLeft.getBand(i);
            float rightAmp = fftRight.getBand(i);
            float amp = (leftAmp + rightAmp) / 2;
            
            if (amp > maxAmp) {
                maxAmp = amp;
                maxBin = i;
            }
        }
        
        if (maxAmp > 0.01) {
            // Convert bin to frequency
            float rawFreq = maxBin * sampleRate / (float)bufferSize;
            float newFreq = constrain(rawFreq, yMin, yMax);
            
            // Smooth the transition
            currentFrequency = prevFrequency + (newFreq - prevFrequency) * FREQ_SMOOTH;
            prevFrequency = currentFrequency;
        }
    }
    
    private float calculateStereoWidth() {
        float totalEnergy = 0;
        float stereoDifference = 0;
        float significantBands = 0;
        float noiseThreshold = 0.05f;
        
        for (int i = 0; i < fftLeft.specSize(); i++) {
            float leftAmp = fftLeft.getBand(i);
            float rightAmp = fftRight.getBand(i);
            float avgAmp = (leftAmp + rightAmp) / 2;
            
            if (avgAmp > noiseThreshold) {
                float difference = abs(leftAmp - rightAmp);
                
                if (avgAmp > 0) {
                    float normalizedDiff = difference / avgAmp;
                    stereoDifference += normalizedDiff * avgAmp;
                    totalEnergy += avgAmp;
                    significantBands++;
                }
            }
        }
        
        if (significantBands > 0 && totalEnergy > noiseThreshold) {
            float newStereoWidth = (stereoDifference / totalEnergy);
            newStereoWidth = prevStereoWidth + (newStereoWidth - prevStereoWidth) * 0.3f;
            prevStereoWidth = newStereoWidth;
            return constrain(newStereoWidth, 0, 1);
        } else {
            prevStereoWidth *= 0.95f;
            return prevStereoWidth;
        }
    }
}

class AudioFileAnalyzer extends BaseAudioAnalyzer {
    private AudioPlayer audioPlayer;
    
    public AudioFileAnalyzer(PApplet parent, float yMin, float yMax, int sampleRate, int bufferSize) {
        super(parent, yMin, yMax, sampleRate, bufferSize);
    }
    
    public void loadFile(String filename) {
        if (audioPlayer != null) {
            audioPlayer.close();
        }
        audioPlayer = minim.loadFile(filename, BUFFER_SIZE);
    }
    
    @Override
    public void play(boolean loop) {
        if (audioPlayer != null) {
            if (loop) {
                audioPlayer.loop();
            } else {
                audioPlayer.play();
            }
        }
    }
    
    @Override
    public void pause() {
        if (audioPlayer != null) {
            audioPlayer.pause();
        }
    }
    
    @Override
    public void update() {
        if (audioPlayer != null && audioPlayer.isPlaying()) {
            bufferLeft = audioPlayer.left.toArray();
            bufferRight = audioPlayer.right.toArray();
            processFFT(bufferLeft, bufferRight);
        }
    }
}

class TestSignalAnalyzer extends BaseAudioAnalyzer implements AudioSignal  {
    private Minim minim;
    private AudioOutput out;
    private FrequencySweeper sweeper;
    private float phase = 0;
    private boolean isPlaying = false;
    
    public TestSignalAnalyzer(PApplet parent, float yMin, float yMax, int sampleRate, int bufferSize) {
        super(parent, yMin, yMax, sampleRate, bufferSize);
        minim = new Minim(parent);
        out = minim.getLineOut(Minim.STEREO, bufferSize);
        sweeper = new FrequencySweeper(yMin, yMax, sampleRate);

        out.addSignal(this);

        play(true);
    }
    
    @Override
    public void play(boolean loop) {
        isPlaying = true;
        sweeper.play();
    }
    
    @Override
    public void pause() {
        isPlaying = false;
        sweeper.pause();
    }

    @Override
    public void update() {
        if (!isPlaying) {
            return;
        }
        
        sweeper.update(frameRate);
        sweeper.generateAudioBuffer(bufferLeft, bufferRight);
        
        // Send audio to output using mix.set()
        //for (int i = 0; i < bufferSize; i++) {
            //out.mix.set(i, (bufferLeft[i] + bufferRight[i]) * 0.5f);
        //}
        
        processFFT(bufferLeft, bufferRight);
    }

    // AudioSignal interface methods
    @Override
    public void generate(float[] buffer) {
        if (!isPlaying) {
            Arrays.fill(buffer, 0);
            return;
        }

        // Sum to mono and apply gain reduction
        float gain = 0.5f;
        for (int i = 0; i < buffer.length; i++) {
            // Average left and right channels
            buffer[i] = (bufferLeft[i] + bufferRight[i]) * 0.5f * gain;
        }
    }

    // AudioSignal interface methods
    @Override
    public void generate(float[] left, float[] right) {
        if (!isPlaying) {
            Arrays.fill(left, 0);
            Arrays.fill(right, 0);
            return;
        }

        arrayCopy(bufferLeft, left);
        arrayCopy(bufferRight, right);
        //left = bufferLeft;
        //right = bufferRight;
        //for (int i = 0; i < left.length; i++) {
        //  left[i] = bufferLeft[i];
        //  right[i] = bufferRight[i];
        //}
    }
}

class FrequencySweeper {
    private float currentFreq;
    private float minFreq;
    private float maxFreq;
    private int sampleRate;
    private float sweepDuration = 1.0;
    private float sweepTime = 0.0;
    private boolean sweepingUp = true;
    private float stereoWidth = 0.0;
    private float phase = 0.0;
    private boolean isPlaying = true;
    
    public FrequencySweeper(float minFreq, float maxFreq, int sampleRate) {
        this.minFreq = minFreq;
        this.maxFreq = maxFreq;
        this.currentFreq = minFreq;
        this.sampleRate = sampleRate;
    }

    public float getCurrentFrequency() { return currentFreq; }
    public float getStereoWidth() { return stereoWidth; }
    public float getPhase() { return phase; }
    public void setPhase(float phase) { this.phase = phase; }

    public void play() {
        isPlaying = true;
    }

    public void pause() {
        isPlaying = false;
    }

    public void update(float frameRate) {
        if(!isPlaying) return;

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

    private void generateAudioBuffer(float[] bufferLeft, float[] bufferRight) {
        float phaseIncrement = TWO_PI * getCurrentFrequency() / sampleRate;
        float phase = getPhase();
        float stereoWidth = getStereoWidth();
        
        for (int i = 0; i < bufferLeft.length; i++) {
            float signal = 0.5 * sin(phase);
            
            bufferLeft[i] = signal * (1.0 - stereoWidth * 0.5);
            bufferRight[i] = signal * (1.0 + stereoWidth * 0.5);
            
            phase += phaseIncrement;
            if (phase > TWO_PI) phase -= TWO_PI;
        }
        
        setPhase(phase);
    }
}

class TimeHistoryManager {
    private TimeSlice[] timeHistory;
    private int currentSliceIndex = 0;
    private float lastUpdateTime = 0;
    private final int TIME_HISTORY = 1000; // milliseconds
    private ParticleGrid grid;
    private float currentStereoWidth = 0.0;
    
    public TimeHistoryManager(Hypercube cube) {
        timeHistory = new TimeSlice[cube.numLayers];
        for (int i = 0; i < cube.numLayers; i++) {
            timeHistory[i] = new TimeSlice(cube.gridSize);
        }
        grid = new ParticleGrid(cube);
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
                float particleSize = Hypercube.POINT_SIZE_MIN;
                if (freqFactor > 0) {
                    particleSize += (Hypercube.POINT_SIZE_MAX - Hypercube.POINT_SIZE_MIN) * freqFactor * stereoFactor;
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

// Helper class for sorting
class ParticlePosition implements Comparable<ParticlePosition> {
    int index;
    float z;
    
    ParticlePosition(int index, float z) {
        this.index = index;
        this.z = z;
    }
    
    @Override
    public int compareTo(ParticlePosition other) {
        return Float.compare(other.z, this.z);
    }
}

class ParticleRenderer {
    private Hypercube cube;
    private ParticleGrid grid;
    private float fieldSeparation = 1.7;
    private PGraphicsOpenGL pgl;
    
    // Color constants
    private static final color TOP_COLOR = #ff0088;    // Light blue
    private static final color BOTTOM_COLOR = #1188ff; // Pink
    
    // Pre-allocated buffers
    private float[] xBuffer;
    private float[] yBuffer;
    private float[] zBuffer;
    private float[] sizeBuffer;
    private float[] alphaBuffer;
    private color[] colorBuffer;
    private int[] sortIndices;
    private int[] staticSortIndices;
    private int particleCount;
    
    // Pre-calculated positions
    private float[] staticXPositions;
    private float[] staticYPositions;
    private float[] staticZPositions;
    
    // Constants
    private static final float OPACITY_THRESHOLD = 0.05f;
    private static final int MAX_PARTICLES = 16384;
    
    public ParticleRenderer(Hypercube cube) {
        this.cube = cube;
        this.grid = new ParticleGrid(cube);
        
        // Pre-allocate buffers
        xBuffer = new float[MAX_PARTICLES];
        yBuffer = new float[MAX_PARTICLES];
        zBuffer = new float[MAX_PARTICLES];
        sizeBuffer = new float[MAX_PARTICLES];
        alphaBuffer = new float[MAX_PARTICLES];
        colorBuffer = new color[MAX_PARTICLES];
        sortIndices = new int[MAX_PARTICLES];
        
        // Initialize OpenGL stuff
        //pgl = (PGraphicsOpenGL)g;
        
        // Pre-calculate positions
        initializeStaticPositions();
        // Pre-calculate sort order
        preCalculateSortOrder();
    }
    
    private color interpolateColor(float t) {
        float r1 = red(TOP_COLOR);
        float g1 = green(TOP_COLOR);
        float b1 = blue(TOP_COLOR);
        
        float r2 = red(BOTTOM_COLOR);
        float g2 = green(BOTTOM_COLOR);
        float b2 = blue(BOTTOM_COLOR);
        
        float r = lerp(r1, r2, t);
        float g = lerp(g1, g2, t);
        float b = lerp(b1, b2, t);
        
        return color(r, g, b);
    }
    
    public void setFieldSeparation(float separation) {
        this.fieldSeparation = separation;
    }

    private void initializeStaticPositions() {
        staticXPositions = new float[cube.gridSize];
        staticYPositions = new float[cube.gridSize];
        staticZPositions = new float[cube.numLayers];
        
        float halfSize = cube.cubeSize / 2.0f;
        
        for (int i = 0; i < cube.gridSize; i++) {
            staticXPositions[i] = map(i, 0, cube.gridSize - 1, -halfSize, halfSize);
            staticYPositions[i] = map(i, 0, cube.gridSize - 1, halfSize, -halfSize);
        }
        
        // Invert the z-positions so newer layers (smaller indices) are closer to camera
        for (int i = 0; i < cube.numLayers; i++) {
            staticZPositions[i] = map(i, 0, cube.numLayers - 1, -halfSize, halfSize);
        }
    }

    private void preCalculateSortOrder() {
        int totalParticles = cube.numLayers * cube.gridSize * cube.gridSize;
        staticSortIndices = new int[totalParticles];
        
        ArrayList<ParticlePosition> positions = new ArrayList<ParticlePosition>();
        
        for (int layer = 0; layer < cube.numLayers; layer++) {
            float z = staticZPositions[layer];
            
            for (int row = 0; row < cube.gridSize; row++) {
                for (int col = 0; col < cube.gridSize; col++) {
                    int index = (layer * cube.gridSize * cube.gridSize) + 
                              (row * cube.gridSize) + col;
                    positions.add(new ParticlePosition(index, z));
                }
            }
        }
        
        // Changed sort order to match new z-axis direction
        Collections.sort(positions, Collections.reverseOrder());
        
        for (int i = 0; i < positions.size(); i++) {
            staticSortIndices[i] = positions.get(i).index;
        }
    }

    private float calculateStereoFieldOpacity(float position, float stereoWidth) {
        float stereoPosition = (position * 2.0f) - 1.0f;
        float peakPosition = stereoPosition < 0 ? -stereoWidth : stereoWidth;
        float distance = abs(stereoPosition - peakPosition);
        return max(0, 1.0f - (distance * fieldSeparation));
    }
    
    private void gatherVisibleParticles(TimeHistoryManager timeHistory, float stereoWidth) {
        particleCount = 0;
        float halfSize = cube.cubeSize / 2.0f;
        
        for (int sortedIdx : staticSortIndices) {
            int totalPerLayer = cube.gridSize * cube.gridSize;
            int layer = sortedIdx / totalPerLayer;
            int remainder = sortedIdx % totalPerLayer;
            int row = remainder / cube.gridSize;
            int col = remainder % cube.gridSize;
            
            TimeSlice slice = timeHistory.getSlice(layer);
            float size = slice.getSizeAt(row, col);
            
            if (size > Hypercube.POINT_SIZE_MIN) {
                float xPos = (float)col / (cube.gridSize - 1);
                float opacity = calculateStereoFieldOpacity(xPos, stereoWidth);
                
                if (opacity > OPACITY_THRESHOLD && particleCount < MAX_PARTICLES) {
                    xBuffer[particleCount] = staticXPositions[col];
                    yBuffer[particleCount] = staticYPositions[row];
                    zBuffer[particleCount] = staticZPositions[layer];
                    sizeBuffer[particleCount] = size;
                    alphaBuffer[particleCount] = opacity * 255;
                    
                    // Calculate color based on y-position
                    float yNormalized = map(staticYPositions[row], -halfSize, halfSize, 0, 1);
                    colorBuffer[particleCount] = interpolateColor(yNormalized);
                    
                    sortIndices[particleCount] = particleCount;
                    particleCount++;
                }
            }
        }
    }
    
    private void renderParticleBatch() {
        if (particleCount == 0) return;

        noStroke();
        
        float lastSize = -1;
        color lastColor = -1;
        float lastAlpha = -1;
        
        beginShape(POINTS);
        for (int i = 0; i < particleCount; i++) {
            int idx = sortIndices[i];
            
            // Only update size and color if they've changed
            float currentSize = sizeBuffer[idx] * 3;
            color currentColor = colorBuffer[idx];
            float currentAlpha = alphaBuffer[idx] * 0.8f;
            
            if (currentSize != lastSize) {
                strokeWeight(currentSize);
                lastSize = currentSize;
            }
            
            if (currentColor != lastColor || currentAlpha != lastAlpha) {
                stroke(red(currentColor), green(currentColor), blue(currentColor), currentAlpha);
                lastColor = currentColor;
                lastAlpha = currentAlpha;
            }
            
            vertex(xBuffer[idx], yBuffer[idx], zBuffer[idx]);
        }
        endShape();
        
        resetShader();
    }
    
    public void render(TimeHistoryManager timeHistory, float stereoWidth) {
        pushStyle();
        
        hint(ENABLE_DEPTH_TEST);
        hint(DISABLE_DEPTH_SORT);
        blendMode(ADD);
        ((PGraphicsOpenGL)g).smooth(4);
        
        gatherVisibleParticles(timeHistory, stereoWidth);
        renderParticleBatch();
        
        blendMode(BLEND);
        popStyle();
    }
}

class ParticleGrid {    
    private float[] xPositions;
    private float[] yPositions;
    private float[] zPositions;
    private float[] rowFrequencies;
    private float[] stereoPositions;
    private float[] centerDistances;
    
    public float getXPosition(int col) { return xPositions[col]; }
    public float getYPosition(int row) { return yPositions[row]; }
    public float getZPosition(int layer) { return zPositions[layer]; }
    public float getRowFrequency(int row) { return rowFrequencies[row]; }
    public float getCenterDistance(int col) { return centerDistances[col]; }

    public ParticleGrid(Hypercube cube) {
        initializeGridPositions(cube);
    }

    private void initializeGridPositions(Hypercube cube) {
        xPositions = new float[cube.gridSize];
        yPositions = new float[cube.gridSize];
        zPositions = new float[cube.numLayers];
        rowFrequencies = new float[cube.gridSize];
        stereoPositions = new float[cube.gridSize];
        centerDistances = new float[cube.gridSize];
        
        float halfSize = cube.cubeSize / 2;
        
        for (int i = 0; i < cube.gridSize; i++) {
            xPositions[i] = map(i, 0, cube.gridSize - 1, -halfSize, halfSize);
            yPositions[i] = map(i, 0, cube.gridSize - 1, halfSize, -halfSize);
            rowFrequencies[i] = cube.yMin * pow(cube.yMax/cube.yMin, (float)i/(cube.gridSize-1));
            stereoPositions[i] = (float)i / (cube.gridSize - 1);
            centerDistances[i] = abs(stereoPositions[i] - 0.5) * 2.0;
        }
        
        // Changed mapping to make particles move away (positive z)
        for (int i = 0; i < cube.numLayers; i++) {
            zPositions[i] = map(i, 0, cube.numLayers - 1, halfSize, -halfSize);
        }
    }
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
