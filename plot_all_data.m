function plot_all_data(data)

% Input validation
required_fields = {'x', 'y', 'z', 'actual_prns', 'satPeaks_before', 'acqThreshold_before', ...
    'satPeaks_after', 'acqThreshold_after', 'X', 'Y', 'Gain_dB', 'doa_el_actual', ...
    'doa_az_actual', 'doa_el_spoof', 'doa_az_spoof', 'u_sp', 'u_jm', 'u_y', 'u_h'};
for i = 1:length(required_fields)
    if ~isfield(data, required_fields{i})
        error('plot_all_data:MissingField', 'Required field "%s" is missing in input data structure.', required_fields{i});
    end
end

% Validate sizes
n = length(data.actual_prns);
if length(data.x) < n+2 || length(data.y) < n+2 || length(data.z) < n+2
    error('plot_all_data:InvalidSize', 'x, y, z must have length at least %d (n+2 for %d satellites plus spoofer and jammer).', n+2, n);
end
if length(data.satPeaks_before) < 32 || length(data.satPeaks_after) < 32
    error('plot_all_data:InvalidSize', 'satPeaks_before and satPeaks_after must have length at least 32.');
end
if ~isscalar(data.acqThreshold_before) || ~isscalar(data.acqThreshold_after)
    error('plot_all_data:InvalidSize', 'acqThreshold_before and acqThreshold_after must be scalars.');
end
if length(data.doa_el_actual) ~= n || length(data.doa_az_actual) ~= n
    error('plot_all_data:InvalidSize', 'doa_el_actual and doa_az_actual must have length %d (same as actual_prns).', n);
end
if ~isscalar(data.doa_el_spoof) || ~isscalar(data.doa_az_spoof)
    error('plot_all_data:InvalidSize', 'doa_el_spoof and doa_az_spoof must be scalars.');
end

% Declare persistent variables for figure and axes handles
persistent fig_handle ax_handles

% Check if the figure exists and is valid; if not, create it
if isempty(fig_handle) || ~ishandle(fig_handle)
    fig_handle = figure('Position', [100 100 1600 900], 'Color', 'w', 'Name', 'All Data Plots', 'NumberTitle', 'off');
    ax_handles = gobjects(2, 3); % Initialize a 2x3 array for axes handles
    ax_handles(1,1) = subplot(2,3,1);
    ax_handles(1,2) = subplot(2,3,2);
    ax_handles(1,3) = subplot(2,3,3);
    ax_handles(2,1) = subplot(2,3,4);
    ax_handles(2,2) = subplot(2,3,5);
    ax_handles(2,3) = subplot(2,3,6);
end

% Bring the figure to the front
figure(fig_handle);

% Subplot 1: 3D Visualization with Satellites, Spoofer, Jammer, and Antenna
axes(ax_handles(1,1));
cla; % Clear the axes
hold on;
plot3(0, 0, 0, 'ro', 'MarkerFaceColor', 'r');
text(0, 0, 0, ' Antenna', 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'center');
quiver3(data.x(1:length(data.actual_prns)), data.y(1:length(data.actual_prns)), data.z(1:length(data.actual_prns)), ...
    -data.x(1:length(data.actual_prns)), -data.y(1:length(data.actual_prns)), -data.z(1:length(data.actual_prns)), 0.5, ...
    'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', 'b');
quiver3(data.x(end-1), data.y(end-1), data.z(end-1), -data.x(end-1), -data.y(end-1), -data.z(end-1), 0.5, ...
    'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', 'm');
quiver3(data.x(end), data.y(end), data.z(end), -data.x(end), -data.y(end), -data.z(end), 0.5, ...
    'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', 'c');
satelliteSize = 200;
for i = 1:length(data.actual_prns)
    scatter3(data.x(i), data.y(i), data.z(i), satelliteSize, 'p', 'filled', 'MarkerFaceColor', 'b');
    text(data.x(i), data.y(i), data.z(i), sprintf('Sat %d', i), 'Color', 'g', 'FontSize', 10, 'HorizontalAlignment', 'center');
end
scatter3(data.x(end-1), data.y(end-1), data.z(end-1), satelliteSize, 'p', 'filled', 'MarkerFaceColor', 'r');
text(data.x(end-1), data.y(end-1), data.z(end-1), 'Spoofer', 'Color', 'r', 'FontSize', 10, 'HorizontalAlignment', 'center');
scatter3(data.x(end), data.y(end), data.z(end), satelliteSize, 'p', 'filled', 'MarkerFaceColor', 'c');
text(data.x(end), data.y(end), data.z(end), 'Jammer', 'Color', 'c', 'FontSize', 10, 'HorizontalAlignment', 'center');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Visualization with Satellites, Spoofer, Jammer, and Antenna');
grid on; axis equal; view(30, 30);
hold off;

% Subplot 2: Histogram for Acquisition Results (Before)
axes(ax_handles(1,2));
cla;
hold on;
AcqnsMetric_before = data.satPeaks_before(:); % Ensure column vector
for i = 1:32
    amplitude = AcqnsMetric_before(i);
    maxAmplitude = max(amplitude);
    if maxAmplitude > data.acqThreshold_before
        bar(i, maxAmplitude, 'FaceColor', 'green');
    else
        bar(i, maxAmplitude, 'FaceColor', 'blue');
    end
end
title('Acquisition Results (Before AS)');
xlabel('Satellite PRN Numbers');
ylabel('Metric Value');
h1 = bar(nan, 'FaceColor', 'blue'); % Dummy bar for 'Not Acquired'
h2 = bar(nan, 'FaceColor', 'green'); % Dummy bar for 'Acquired'
legend([h1 h2], {'Not Acquired (Before AS)', 'Acquired (Before AS)'});
hold off;

% Subplot 3: Histogram for Acquisition Results (After)
axes(ax_handles(1,3));
cla;
hold on;
AcqnsMetric_after = data.satPeaks_after(:); % Ensure column vector
for i = 1:32
    amplitude = AcqnsMetric_after(i);
    maxAmplitude = max(amplitude);
    if maxAmplitude > data.acqThreshold_after
        bar(i, maxAmplitude, 'FaceColor', 'green');
    else
        bar(i, maxAmplitude, 'FaceColor', 'blue');
    end
end
title('Acquisition Results (After AS)');
xlabel('Satellite PRN Numbers');
ylabel('Metric Value');
h1 = bar(nan, 'FaceColor', 'blue'); % Dummy bar for 'Not Acquired'
h2 = bar(nan, 'FaceColor', 'green'); % Dummy bar for 'Acquired'
legend([h1 h2], {'Not Acquired (After AS)', 'Acquired (After AS)'});
hold off;

% Subplot 4: Polar Beam Pattern
axes(ax_handles(2,1));
cla;
p = pcolor(data.X, data.Y, data.Gain_dB);
set(p, 'EdgeColor', 'none');
colormap(jet);
caxis([-60 0]);
axis equal off;
hold on;
for elc = [0 30 60]
    rc = 90 - elc;
    th = linspace(0, 2*pi, 360);
    plot(rc * cos(th), rc * sin(th), 'k:', 'LineWidth', 1);
    text((rc + 2) * cosd(135), (rc + 2) * sind(135), sprintf('%d° el', elc), ...
         'FontSize', 9, 'Color', 'k');
end
for az = 0:30:330
    r0 = 92;
    text(r0 * cosd(az), r0 * sind(az), sprintf('%d°', az), ...
         'FontSize', 9, 'HorizontalAlignment', 'center');
end
for i = 1:length(data.actual_prns)
    ri = 90 - data.doa_el_actual(i);
    xi = ri * cosd(data.doa_az_actual(i));
    yi = ri * sind(data.doa_az_actual(i));
    plot(xi, yi, 'gp', 'MarkerSize', 3, 'MarkerFaceColor', 'g');
    text(xi, yi, sprintf('%d_A', data.actual_prns(i)), 'Color', 'g', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
rs = 90 - data.doa_el_spoof;
xs = rs * cosd(data.doa_az_spoof);
ys = rs * sind(data.doa_az_spoof);
plot(xs, ys, 'rp', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
text(xs, ys, 'S', 'Color', 'b', 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
hcb = colorbar('SouthOutside');
hcb.Label.String = 'Relative Gain (dB)';
title('Polar Beam Pattern (Null-Steered) vs Azimuth & Elevation', 'FontSize', 12);
hold off;

% Subplot 5: 3D DOA Visualization
axes(ax_handles(2,2));
cla;
hold on;
axis equal; grid on; view(30,30);
plot3(0,0,0,'ko','MarkerFaceColor','y');
text(0,0,0,'Antenna','HorizontalAlignment','center');
L = 0.6;
quiver3(0,0,0, -data.u_sp(1)*L, -data.u_sp(2)*L, -data.u_sp(3)*L, 'm','LineWidth',2);
text(-data.u_sp(1)*L, -data.u_sp(2)*L, -data.u_sp(3)*L, 'Orig Spoofer','Color','m');
quiver3(0,0,0, -data.u_jm(1)*L, -data.u_jm(2)*L, -data.u_jm(3)*L, 'c','LineWidth',2);
text(-data.u_jm(1)*L, -data.u_jm(2)*L, -data.u_jm(3)*L, 'Orig Jammer','Color','c');
quiver3(0,0,0, -data.u_y(1)*L, -data.u_y(2)*L, -data.u_y(3)*L, 'r--','LineWidth',2);
text(-data.u_y(1)*L, -data.u_y(2)*L, -data.u_y(3)*L, 'Estimated y','Color','r');
quiver3(0,0,0, -data.u_h(1)*L, -data.u_h(2)*L, -data.u_h(3)*L, 'g-.','LineWidth',2);
text(-data.u_h(1)*L, -data.u_h(2)*L, -data.u_h(3)*L, 'h vector','Color','g');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D DOA: Spoofer, Jammer, Estimated y & h');
hold off;



% Subplot 6: after pmax 
axes(ax_handles(2,3));
cla;
 hold on;
AcqnsMetric_after_p = data.satPeaks_after_p(:); % Ensure column vector
for i = 1:32
    amplitude = AcqnsMetric_after_p(i);
    maxAmplitude = max(amplitude);
    if maxAmplitude > data.acqThreshold_after_p
        bar(i, maxAmplitude, 'FaceColor', 'green');
    else
        bar(i, maxAmplitude, 'FaceColor', 'blue');
    end
end
title('Acquisition Results (After AS with Pmax)');
xlabel('Satellite PRN Numbers');
ylabel('Metric Value');
h1 = bar(nan, 'FaceColor', 'blue'); % Dummy bar for 'Not Acquired'
h2 = bar(nan, 'FaceColor', 'green'); % Dummy bar for 'Acquired'
legend([h1 h2], {'Not Acquired (After AS with Pmax)', 'Acquired (After AS with Pmax)'});
% hold off;

end